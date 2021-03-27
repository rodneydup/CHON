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
  Spring(const char* name) {
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
