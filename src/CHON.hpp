#include "./Roboto-Medium.hpp"
#include "./structs.hpp"
#include "Gamma/DFT.h"
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
using namespace al;

class CHON : public App {
 public:
  /**
   * @brief Initilialize the synth interface.
   */
  virtual void onInit() override;

  /**
   * @brief Run once on starup.
   */
  virtual void onCreate() override;

  /**
   * @brief Audio rate processing of synth.
   */
  virtual void onSound(al::AudioIOData &io) override;

  void drawAudioIO(al::AudioIO *io);

  virtual void onAnimate(double dt) override;

  /**
   * @brief Draw rate processing of synth interface.
   */
  virtual void onDraw(al::Graphics &g) override;

  virtual void onResize(int width, int height) override;

  virtual Vec3d unproject(Vec3d screenPos);

  virtual Rayd getPickRay(int screenX, int screenY);

  virtual bool onMouseMove(const Mouse &m) override;

  virtual bool onMouseDown(const al::Mouse &m) override;

  virtual bool onMouseDrag(const Mouse &m) override;

  virtual bool onMouseUp(const Mouse &m) override;

  virtual bool onKeyDown(al::Keyboard const &k) override;

  virtual bool onKeyUp(al::Keyboard const &k) override;

  virtual void chonReset();

  /** MIDI Stuff **/
  //   void initMIDI();
  //   void updateActiveMIDIParams(const al::MIDIMessage &m);

  /**
   * @brief Called everytime a MIDI message is sent.
   */
  //   virtual void onMIDIMessage(const al::MIDIMessage &m) override;

  //   virtual void onExit() override;

  int xParticles = 4;     // Parameter for x particles count
  int yParticles = 1;     // Parameter for y particles count
  al::Mesh particleMesh;  // mesh for drawing particle

  Particle *rightClickedParticle;
  bool isRightClickedParticle = false;
  float rightClickedFreq = 0;
  float rightClickedAmplitude = 1;
  int rightClickedStep = 1;

  ParticleNetwork particleNetwork{xParticles, yParticles};  // 2D vector containing our Particles

  Texture texBlur;  // blurring filter for graph
  int amCounter;    // AM incrementor
  std::mutex resetLock;
  float w = 2;
  float h = 2;

  Reverb<float> reverb;
  gam::NoiseWhite<> tick;
  bool drawGUI = 1;

  /*
  GUI variables
  */
  std::unique_ptr<ChonBundle> xSpringGUI;
  std::unique_ptr<ChonBundle> ySpringGUI;
  Parameter mass{"Mass", "physics", 1.0f, "", 1.0f, 100.0f};      // Master mass
  Parameter damping{"Damping", "physics", 0.0f, "", 0.0f, 3.0f};  // damping
  ParameterBool xFree{"X axis", "Degrees of Freedom", 1};
  ParameterBool yFree{"Y axis", "Degrees of Freedom", 0};
  ParameterBool zFree{"Z axis", "Degrees of Freedom", 0};
  ParameterBool pause{"Pause (press p)", "physics", 0};  // Pause Simulation

  ParameterBool DrawGraph{"Draw Graph", "draw", 1};               // Toggle Drawing Graph
  Parameter graphSpread{"Spread", "draw", 0.0f, "", 0.0f, 1.0f};  // Graph spread
  Parameter graphSpeed{"Speed", "draw", 10.0f, "", 1.0f, 30.0f};  // Graph draw speed
  ParameterMenu graphAxis{"Graph Axis"};  // choose which axis displacement to draw
  ParameterBool drawParticles{"Draw Particles", "draw", 1};    // Toggle drawing particles
  ParameterBool drawBoundaries{"Draw Boundaries", "draw", 0};  // Toggle drawing boundaries

  ParameterBool stereoOn{"Stereo", "Synthesis", 0};                     // Stereo Mode toggle
  ParameterMenu tuningScale{"Scale"};                                   // Tuning of Bell synth
  Parameter tuningRoot{"Root", "Synthesis", 60.0f, "", 1.0f, 1000.0f};  // Root of bell tuning
  ParameterBool additiveSynthOn{"Additive Synth On", "Synthesis", 0};   // Additive Synth toggle
  ParameterBool bellSynthOn{"Bell Synth On", "Synthesis", 0};           // Bell Synth toggle
  ParameterMenu bellAxis{"Axis##bell"};

  Parameter bellVolume{"Volume##bell", "Synthesis", 0.5f, "", 0.0f, 1.0f};  // Volume of bell synth

  Parameter additiveVolume{
    "Volume##additive", "Synthesis", 0.5f, "", 0.0f, 1.0f};  // Volume of bell synth

  ParameterBool fm{"FM On", "Synthesis", 0};  // FM toggle
  ParameterMenu fmAxis{"Axis##FM"};
  Parameter fmFreq{
    "Frequency##FM", "Synthesis", 1.5f, "", 0.0f, 2.0f};              // FM freq (ratio to carrier)
  Parameter fmWidth{"Width##FM", "Synthesis", 2.0f, "", 0.1f, 5.0f};  // FM Width

  ParameterBool am{"AM On", "Synthesis", 0};  // AM toggle
  ParameterMenu amAxis{"Axis##AM"};

  ParameterBool reverbOn{"Reverb On", "Synthesis", 0};               // Reverb
  Parameter reverbTail{"Tail", "Synthesis", 0.15f, "", 0.0f, 1.0f};  // Reverb decay time

  ParameterBool inputOn{"Input On", "Audio", 0};
  ParameterMenu inputMode{"Input Mode"};
  ParameterBool driveStereoSplit{"Stereo Split", "Audio", 0};
  ParameterInt driveParticleXLeft{"X##Left", "Audio", 1, 1, 2};
  ParameterInt driveParticleYLeft{"Y##Left", "Audio", 1, 1, 2};
  ParameterMenu driveAxisLeft{"Drive Axis##Left"};
  ParameterInt driveParticleXRight{"X##Right", "Audio", 1, 1, 2};
  ParameterInt driveParticleYRight{"Y##Right", "Audio", 1, 1, 2};
  ParameterMenu driveAxisRight{"Drive Axis##Right"};
  Parameter inputThreshold{"Input Threshold", "Audio", 1.0f, 0.0f, 2.0f};
  Parameter inputScale{"Input Scaling", "Audio", 1.0f, 0.1f, 2.0f};
  ParameterInt rmsSize{"RMS Samples", "Audio", 2048, 512, 4096};

  /*
  Scales
  */
  std::vector<float> majScale{1.000000, 1.125000, 1.250000, 1.333333,
                              1.500000, 1.666667, 1.875000, 2.000000};
  std::vector<float> pentScale{1.000000, 1.125000, 1.250000, 1.500000, 1.666667, 2.000000};
  std::vector<float> chromScale{1.000000, 1.066667, 1.125000, 1.285714, 1.250000,
                                1.333333, 1.406250, 1.500000, 1.600000, 1.666667,
                                1.750000, 1.875000, 2.000000};
  std::vector<float> bpScale{1.000000, 1.080000, 1.190476, 1.285714, 1.400000, 1.530612, 1.666667,
                             1.800000, 1.960000, 2.142857, 2.333333, 2.520000, 2.777778, 3.000000};

  /*
  OSC
  */
  int port = 16447;             // osc port
  char addr[10] = "127.0.0.1";  // ip address
  osc::Send client;             // create an osc client
  ParameterBool oscOn{"OSC On", "OSC", 0};
  ParameterBool oscX{"X disp", "OSC", 0};
  ParameterBool oscY{"Y disp", "OSC", 0};
  ParameterBool oscZ{"Z disp", "OSC", 0};
  ParameterBool oscPan{"Pan", "OSC", 0};

  void resetOSC() {
    client.open(port, addr);
    std::cout << "New OSC port Selected: \n" + port << std::endl;
    std::cout << addr << std::endl;
  }

  ImGuiWindowFlags flags = ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoMove |
                           ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings |
                           ImGuiWindowFlags_AlwaysAutoResize | ImGuiCond_Once;
  ImFont *bodyFont;
  ImFont *titleFont;

  /*
  Audio IO
  */
  std::string currentAudioDeviceOut;
  std::string currentAudioDeviceIn;

  static const int MAX_AUDIO_OUTS = 2;
  static const int MAX_AUDIO_INS = 2;

  std::array<unsigned int, MAX_AUDIO_OUTS> AudioChanIndexOut;
  std::array<unsigned int, MAX_AUDIO_INS> AudioChanIndexIn;

  bool isPaused = false;
  const int SAMPLE_RATE = 48000;
  double globalSamplingRate = SAMPLE_RATE;
  const int BLOCK_SIZE = 1024;

  int getLeadChannelOut() const { return AudioChanIndexOut[0]; }
  int getLeadChannelIn() const { return AudioChanIndexIn[0]; }

  void setOutChannels(int lead_channel, int max_possible_channels) {
    AudioChanIndexOut[0] = lead_channel;
    if (max_possible_channels == 1) {
      for (int i = 1; i < MAX_AUDIO_OUTS; i++) {
        AudioChanIndexOut[i] = lead_channel;
      }
    } else {
      // assert(lead_channel + (consts::MAX_AUDIO_OUTS) < max_possible_channels);
      for (int i = 1; i < MAX_AUDIO_OUTS; i++) {
        AudioChanIndexOut[i] = lead_channel + i;
      }
    }
  }
  void setInChannels(int lead_channel, int max_possible_channels) {
    AudioChanIndexIn[0] = lead_channel;
    if (max_possible_channels == 1) {
      for (int i = 1; i < MAX_AUDIO_OUTS; i++) {
        AudioChanIndexIn[i] = lead_channel;
      }
    } else {
      // assert(lead_channel + (consts::MAX_AUDIO_OUTS) < max_possible_channels);
      for (int i = 1; i < MAX_AUDIO_OUTS; i++) {
        AudioChanIndexIn[i] = lead_channel + i;
      }
    }
  }
  int getSampleRateIndex() {
    unsigned s_r = (unsigned)globalSamplingRate;
    switch (s_r) {
      case 44100:
        return 0;
      case 48000:
        return 1;
      case 88200:
        return 2;
      case 96000:
        return 3;
      default:
        return 0;
    }
  }

  // Audio input buffers
  unsigned int inBufferSize = 4096;
  RingBuffer inLeft{inBufferSize};
  RingBuffer inRight{inBufferSize};
  float driveForceLeft = 0;
  float driveForceRight = 0;
  gam::STFT stft{
    256,          // Window size
    256,          // Hop size; number of samples between transforms
    0,            // Pad size; number of zero-valued samples appended to window
    gam::HANN,    // Window type: BARTLETT, BLACKMAN, BLACKMAN_HARRIS,
                  //		HAMMING, HANN, WELCH, NYQUIST, or RECTANGLE
    gam::COMPLEX  // Format of frequency samples:
                  //		COMPLEX, MAG_PHASE, or MAG_FREQ
  };
  std::array<float, 129> fftBuffer = {};
  float fftDivision = 1;
  int fftIterator = 0;
};
