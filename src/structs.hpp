#include <array>

#include "Gamma/Filter.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
#include "al/app/al_App.hpp"
#include "al/graphics/al_Shapes.hpp"
#include "al/math/al_Ray.hpp"
#include "al/sound/al_Reverb.hpp"
#include "al/ui/al_BoundingBox.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al/ui/al_Pickable.hpp"

template <class T>  // The class is templated to allow a variety of data types
class SmoothValue {
 public:
  SmoothValue(float initialTime =
                50.0f) {  // how much time, in ms, does it take to arrive (or approach target value)
    arrivalTime = initialTime;  // store the input value (default to 50ms)
    calculateCoefficients();    // calculate a and b (it is a lowpass filter)
    z = 0.0f;
  }

  inline T process() {  // this function will be called per-sample on the data
    z = (targetValue * b) +
        (z * a);  // when evaluated, z is the 'previous z' (it is leftover from last execution)
    return z;     // return the new z value (the output sample)
  }
  void setTime(float newTime) {  // set time (in ms)
    arrivalTime = newTime;       // store the input value
    calculateCoefficients();     // calculate a and b
  }
  void setTarget(T newTarget) { targetValue = newTarget; }
  T getCurrentValue() { return z; }
  T getTargetValue() { return targetValue; }
  float getTime() { return arrivalTime; }

 private:
  void calculateCoefficients() {  // called only when 'setTime' is called (and in constructor)
    a = std::exp(-(M_PI * 2) /
                 (arrivalTime * 0.001f * gam::sampleRate()));  // rearranged lpf coeff calculations
    b = 1.0f - a;
  }
  T targetValue;      // what is the destination (of type T, determind by implementation)
  T currentValue;     // how close to the destination? (output value)
  float arrivalTime;  // how long to take
  float a, b;         // coefficients
  T z;                // storage for previous value
};

struct Spring {
  Spring(const char *name) {
    bundle.name(name);
    bundle << k;
  }
  al::Parameter k{"Stiffness", "physics", 1.0f, "", 1.0f, 50.0f};
  al::ParameterBundle bundle;
};

struct Particle {
 public:
  double velocity[3] = {0, 0, 0};
  double acceleration[3] = {0, 0, 0};
  double equilibrium[3] = {0, 0, 0};
  double displacement[3] = {0, 0, 0};
  double prevDisplacement[3] = {0, 0, 0};
  std::string oscName[3] = {"", "", ""};
  SmoothValue<float> amSmooth;
  SmoothValue<float> fmSmooth;

  float mass;
  al::PickableBB particle;
  al::Mesh graph;
  gam::Sine<> oscillator;
  gam::Sine<> FM;
  gam::Sine<> bell;
  gam::Biquad<> Lop;
  double bellEnv;
  bool zeroTrigger[3] = {0, 0, 0};
  Particle() {
    Lop.type(gam::FilterType(gam::LOW_PASS));
    Lop.set(400, 1);
    oscillator.freq(0);
    bell.freq(0);
    amSmooth.setTime(40.0f);
    fmSmooth.setTime(40.0f);
  }

  void setPos(double x, double y, double z) { particle.pose.setPos(al::Vec3d(x, y, z)); }
  double x() { return particle.pose.get().x(); }
  double y() { return particle.pose.get().y(); }
  double z() { return particle.pose.get().z(); }
  void x(double newX) {
    al::Vec3d position = {newX, particle.pose.get().y(), particle.pose.get().z()};
    particle.pose.setPos(position);
  }
  void y(double newY) {
    al::Vec3d position = {particle.pose.get().x(), newY, particle.pose.get().z()};
    particle.pose.setPos(position);
  }
  void z(double newZ) {
    al::Vec3d position = {particle.pose.get().x(), particle.pose.get().y(), newZ};
    particle.pose.setPos(position);
  }
  void updateDisplacement() {
    for (int i = 0; i < 3; i++) {
      prevDisplacement[i] = displacement[i];
      displacement[i] = particle.pose.get().pos()[i] - equilibrium[i];
    }
  }
  void addVelocity() {
    al::Vec3d position = {this->x() + this->velocity[0], this->y() + this->velocity[1],
                          this->z() + this->velocity[2]};
    particle.pose.setPos(position);
  }
  float getBellEnv() { return Lop(this->bellEnv); }
};

class RingBuffer {
 public:
  RingBuffer(unsigned maxSize) : mMaxSize(maxSize) {
    mBuffer.resize(mMaxSize);
    mTail = -1;
    mPrevSample = 0;
  }

  unsigned getMaxSize() const { return mMaxSize; }

  void resize(unsigned maxSize) {
    mMaxSize = maxSize;
    mBuffer.resize(mMaxSize);
  }

  void push_back(float value) {
    mMutex.lock();
    mTail = (mTail + 1) % mMaxSize;
    mBuffer[mTail] = value;
    mMutex.unlock();
  }

  unsigned getTail() const { return mTail; }

  float at(unsigned index) {
    if (index >= mMaxSize) {
      std::cerr << "RingBuffer index out of range." << std::endl;
      index = index % mMaxSize;
    }
    if (mMutex.try_lock()) {
      mPrevSample = mBuffer.at(index);
      mMutex.unlock();
    }
    return mPrevSample;
  }

  float operator[](unsigned index) { return this->at(index); }

  const float *data() { return mBuffer.data(); }

  float getRMS(unsigned lookBackLength) {
    int start = mTail - lookBackLength;
    if (start < 0) start = mMaxSize + start;

    float val = 0.0;
    for (unsigned i = 0; i < lookBackLength; i++) {
      val += pow(mBuffer[(start + i) % mMaxSize], 2);
    }
    return sqrt(val / lookBackLength);
  }

  void print() const {
    for (auto i = mBuffer.begin(); i != mBuffer.end(); ++i) std::cout << *i << " ";
    std::cout << "\n";
  }

 private:
  std::vector<float> mBuffer;
  unsigned mMaxSize;
  int mTail;
  float mPrevSample;

  std::mutex mMutex;
};

/// BundleGUIManager copied from Allolib and modified to give access to "global"
class ChonBundle {
  /// @ingroup UI
 public:
  ChonBundle(bool global) { mBundleGlobal = global; }

  void drawBundleGUI() {
    std::unique_lock<std::mutex> lk(mBundleLock);
    std::string suffix = "##_bundle_" + mName;
    al::ParameterGUI::drawBundleGroup(mBundles, suffix, mCurrentBundle, mBundleGlobal);
  }

  ChonBundle &registerParameterBundle(al::ParameterBundle &bundle) {
    if (mName.size() == 0 || bundle.name() == mName) {
      std::unique_lock<std::mutex> lk(mBundleLock);
      if (mName.size() == 0) {
        mName = bundle.name();
      }
      mBundles.push_back(&bundle);
    } else {
      std::cout << "Warning: bundle name mismatch. Bundle '" << bundle.name() << "' ingnored."
                << std::endl;
    }
    return *this;
  };

  /// Register parameter using the streaming operator.
  ChonBundle &operator<<(al::ParameterBundle &newBundle) {
    return registerParameterBundle(newBundle);
  }

  /// Register parameter using the streaming operator.
  ChonBundle &operator<<(al::ParameterBundle *newBundle) {
    return registerParameterBundle(*newBundle);
  }

  std::string name() { return mName; }

  int &currentBundle() { return mCurrentBundle; }
  bool &bundleGlobal() { return mBundleGlobal; }
  bool &bundleGlobal(bool global) { mBundleGlobal = global; }
  std::vector<al::ParameterBundle *> bundles() { return mBundles; }

 private:
  std::mutex mBundleLock;
  std::vector<al::ParameterBundle *> mBundles;
  std::string mName;
  int mCurrentBundle{0};
  bool mBundleGlobal{false};
};
