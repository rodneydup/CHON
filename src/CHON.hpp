#include "../external/nativefiledialog/src/include/nfd.h"
#include "./Roboto-Medium.hpp"
#include "./consts.hpp"
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
#include "al_ext/soundfile/al_OutputRecorder.hpp"
#include "nlohmann/json.hpp"
using json = nlohmann::json;
#ifdef _WIN32
#include <stdlib.h>
#define strcasecmp _stricmp
#define PATH_MAX 1024
#include <direct.h>
#define GetCurrentDir _getcwd
#elif __linux__
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#define GetCurrentDir getcwd
#else
#include <limits.h>
#include <mach-o/dyld.h>
#include <stdlib.h>
#define GetCurrentDir getcwd
#endif

#include <stdio.h> /* defines FILENAME_MAX */

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

  /**
   * @brief Modified version of al's soundfilerecordGUI.
   *
   * @param[in] Output recorder object.
   *
   * @param[in] Path of directory where the outputted sound files will be
   * stored.
   *
   * @param[in] Frame rate of outputted file.
   *
   * @param[in] Number of channels of outputted file.
   *
   * @param[in] Amount of space allocated for sound.
   */
  void drawRecorderWidget(al::OutputRecorder *recorder, double frameRate, uint32_t numChannels,
                          std::string directory = "", uint32_t bufferSize = 0);
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

  virtual void onExit() override;

  virtual void chonReset();

  /** MIDI Stuff **/
  //   void initMIDI();
  //   void updateActiveMIDIParams(const al::MIDIMessage &m);

  /**
   * @brief Called everytime a MIDI message is sent.
   */
  //   virtual void onMIDIMessage(const al::MIDIMessage &m) override;

  //   virtual void onExit() override;

  int xParticles = 4;  // Parameter for x particles count
  int yParticles = 1;  // Parameter for y particles count
  int zParticles = 1;  // Parameter for z particles count

  al::Mesh particleMesh;  // mesh for drawing particle

  Particle *rightClickedParticle;
  bool isRightClickedParticle = false;
  float rightClickedFreq = 0;
  float rightClickedAmplitude = 1;
  int rightClickedStep = 1;

  ParticleNetwork particleNetwork{xParticles, yParticles,
                                  zParticles};  // 3D vector containing our Particles

  Texture texBlur;  // blurring filter for graph
  int amCounter;    // AM incrementor
  std::mutex resetLock;
  float w = 2;
  float h = 2;

  Reverb<float> reverb;
  al::OutputRecorder mRecorder;
  std::string soundOutput, execDir, execPath, userPath, configFile, presetsPath;
  nfdresult_t result;
  nfdchar_t *outPath = NULL;
  bool freezedt = false;  // variable needed because NFD causes huge delta time
                          // which messes with simulation
  std::unique_ptr<al::PresetHandler> mPresets;

  float windowWidth, windowHeight;
  bool isFirstLaunch;
  bool saveDefaultAudio = false;

  bool drawGUI = 1;

  /*
  GUI variables
  */
  std::unique_ptr<ChonBundle> xSpringGUI;
  std::unique_ptr<ChonBundle> ySpringGUI;
  std::unique_ptr<ChonBundle> zSpringGUI;
  Parameter mass{"Mass", "physics", 1.0f, 1.0f, 100.0f};      // Master mass
  Parameter damping{"Damping", "physics", 0.0f, 0.0f, 3.0f};  // damping
  ParameterBool xFree{"X axis", "Degrees of Freedom", 1};
  ParameterBool yFree{"Y axis", "Degrees of Freedom", 0};
  ParameterBool zFree{"Z axis", "Degrees of Freedom", 0};
  ParameterBool pause{"Pause (press p)", "physics", 0};  // Pause Simulation

  ParameterBool DrawGraph{"Draw Graph", "draw", 1};           // Toggle Drawing Graph
  Parameter graphSpread{"Spread", "draw", 0.0f, 0.0f, 1.0f};  // Graph spread
  Parameter graphSpeed{"Speed", "draw", 10.0f, 1.0f, 30.0f};  // Graph draw speed
  ParameterMenu graphAxis{"Graph Axis"};  // choose which axis displacement to draw
  ParameterBool drawParticles{"Draw Particles", "draw", 1};    // Toggle drawing particles
  ParameterBool drawBoundaries{"Draw Boundaries", "draw", 0};  // Toggle drawing boundaries

  ParameterBool stereoOn{"Stereo", "Synthesis", 0};                    // Stereo Mode toggle
  ParameterMenu tuningScale{"Scale"};                                  // Tuning of Bell synth
  Parameter tuningRoot{"Root", "Synthesis", 60.0f, 1.0f, 1000.0f};     // Root of bell tuning
  ParameterBool additiveSynthOn{"Additive Synth On", "Synthesis", 0};  // Additive Synth toggle
  ParameterBool bellSynthOn{"Bell Synth On", "Synthesis", 0};          // Bell Synth toggle
  ParameterMenu bellAxis{"Axis##bell"};

  Parameter bellVolume{"Volume##bell", "Synthesis", 0.5f, 0.0f, 1.0f};  // Volume of bell synth

  Parameter additiveVolume{"Volume##additive", "Synthesis", 0.5f, 0.0f,
                           1.0f};  // Volume of bell synth

  ParameterBool fm{"FM On", "Synthesis", 0};  // FM toggle
  ParameterMenu fmAxis{"Axis##FM"};
  Parameter fmFreq{"Frequency##FM", "Synthesis", 1.5f, 0.0f, 2.0f};  // FM freq (ratio to carrier)
  Parameter fmWidth{"Width##FM", "Synthesis", 2.0f, 0.1f, 5.0f};     // FM Width

  ParameterBool am{"AM On", "Synthesis", 0};  // AM toggle
  ParameterMenu amAxis{"Axis##AM"};

  ParameterBool reverbOn{"Reverb On", "Synthesis", 0};           // Reverb
  Parameter reverbTail{"Tail", "Synthesis", 0.15f, 0.0f, 1.0f};  // Reverb decay time

  ParameterBool inputOn{"Input On", "Audio", 0};
  ParameterMenu inputMode{"Input Mode"};
  ParameterBool driveStereoSplit{"Stereo Split", "Audio", 0};
  ParameterInt driveParticleXLeft{"X##Left", "Audio", 1, 1, 2};
  ParameterInt driveParticleYLeft{"Y##Left", "Audio", 1, 1, 2};
  ParameterInt driveParticleZLeft{"Z##Left", "Audio", 1, 1, 2};
  ParameterMenu driveAxisLeft{"Drive Axis##Left"};
  ParameterInt driveParticleXRight{"X##Right", "Audio", 1, 1, 2};
  ParameterInt driveParticleYRight{"Y##Right", "Audio", 1, 1, 2};
  ParameterInt driveParticleZRight{"Y##Right", "Audio", 1, 1, 2};
  ParameterMenu driveAxisRight{"Drive Axis##Right"};
  Parameter inputThreshold{"Input Threshold", "Audio", 1.0f, 0.0f, 2.0f};
  Parameter inputScale{"Input Scaling", "Audio", 1.0f, 0.1f, 10.0f};
  Parameter FFTlog{"Logarithmic Scaling", "Audio", 2.0f, 1.0f, 10.0f};
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
  char addr[15] = "127.0.0.1";  // ip address
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

  std::string buf1 = "test.wav";
  ImGuiInputTextCallback buf1Callback;
  void *buf1CallbackUserData;
  std::string recordFilename;

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

  std::array<unsigned int, consts::MAX_AUDIO_OUTS> AudioChanIndexOut;
  std::array<unsigned int, consts::MAX_AUDIO_INS> AudioChanIndexIn;

  bool isPaused = false;
  double globalSamplingRate = consts::SAMPLE_RATE;
  double globalBufferSize = consts::BUFFER_SIZE;

  int getLeadChannelOut() const { return AudioChanIndexOut[0]; }
  int getLeadChannelIn() const { return AudioChanIndexIn[0]; }

  void setOutChannels(int lead_channel, int max_possible_channels) {
    AudioChanIndexOut[0] = lead_channel;
    if (max_possible_channels == 1) {
      for (int i = 1; i < consts::MAX_AUDIO_OUTS; i++) {
        AudioChanIndexOut[i] = lead_channel;
      }
    } else {
      // assert(lead_channel + (consts::MAX_AUDIO_OUTS) < max_possible_channels);
      for (int i = 1; i < consts::MAX_AUDIO_OUTS; i++) {
        AudioChanIndexOut[i] = lead_channel + i;
      }
    }
  }
  void setInChannels(int lead_channel, int max_possible_channels) {
    AudioChanIndexIn[0] = lead_channel;
    if (max_possible_channels == 1) {
      for (int i = 1; i < consts::MAX_AUDIO_INS; i++) {
        AudioChanIndexIn[i] = lead_channel;
      }
    } else {
      // assert(lead_channel + (consts::MAX_AUDIO_OUTS) < max_possible_channels);
      for (int i = 1; i < consts::MAX_AUDIO_INS; i++) {
        AudioChanIndexIn[i] = lead_channel + i;
      }
    }
  }
  int getSampleRateIndex() {
    unsigned s_r = (unsigned)getGlobalSamplingRate();
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

  int getBufSizeIndex() {
    unsigned bufSize = (int)getGlobalBufferSize();
    switch (bufSize) {
      case 128:
        return 0;
      case 256:
        return 1;
      case 512:
        return 2;
      case 1024:
        return 3;
      case 2048:
        return 4;
      default:
        return 0;
    }
  }

  void setGlobalSamplingRate(float sampling_rate) { globalSamplingRate = sampling_rate; }
  double getGlobalSamplingRate() { return globalSamplingRate; }

  void setGlobalBufferSize(float buffer_size) { globalBufferSize = buffer_size; }
  double getGlobalBufferSize() { return globalBufferSize; }

  void setSoundOutputPath(std::string sound_output_path) {
    soundOutput = al::File::conformPathToOS(sound_output_path);
  }
  void setAudioSettings(float sample_rate, int buffer_size) {
    globalSamplingRate = sample_rate;
    configureAudio(sample_rate, buffer_size, consts::MAX_AUDIO_OUTS, -1);
  }
  void setWindowDimensions(float width, float height) {
    windowWidth = width;
    windowHeight = height;
    dimensions(width, height);
  }
  void setFirstLaunch(bool is_first_launch) { isFirstLaunch = is_first_launch; }
  void setAudioDeviceOut(std::string audio_device) { currentAudioDeviceOut = audio_device; }
  void setAudioDeviceIn(std::string audio_device) { currentAudioDeviceIn = audio_device; }

  // JSON config file stuff
  bool initJsonConfig() {
    json config;
    std::ifstream ifs(userPath + configFile);

    if (ifs.is_open()) {
      config = json::parse(ifs);

      if (config.find(consts::SOUND_OUTPUT_PATH_KEY) == config.end())
        config[consts::SOUND_OUTPUT_PATH_KEY] =
          al::File::conformPathToOS(userPath + consts::DEFAULT_SOUND_OUTPUT_PATH);

      if (config.find(consts::SAMPLE_RATE_KEY) == config.end())
        config[consts::SAMPLE_RATE_KEY] = consts::SAMPLE_RATE;

      if (config.find(consts::BUFFER_SIZE_KEY) == config.end())
        config[consts::BUFFER_SIZE_KEY] = consts::BUFFER_SIZE;

      // if (config.find(con::FONT_SCALE_KEY) == config.end())
      //   config[con::FONT_SCALE_KEY] = con::FONT_SCALE;

      if (config.find(consts::WINDOW_WIDTH_KEY) == config.end())
        config[consts::WINDOW_WIDTH_KEY] = consts::WINDOW_WIDTH;

      if (config.find(consts::WINDOW_HEIGHT_KEY) == config.end())
        config[consts::WINDOW_HEIGHT_KEY] = consts::WINDOW_HEIGHT;

      if (config.find(consts::IS_FIRST_LAUNCH_KEY) == config.end())
        config[consts::IS_FIRST_LAUNCH_KEY] = consts::IS_FIRST_LAUNCH;

      if (config.find(consts::DEFAULT_AUDIO_DEVICE_OUT_KEY) == config.end())
        config[consts::DEFAULT_AUDIO_DEVICE_OUT_KEY] = consts::DEFAULT_AUDIO_DEVICE_OUT;

      if (config.find(consts::DEFAULT_AUDIO_DEVICE_IN_KEY) == config.end())
        config[consts::DEFAULT_AUDIO_DEVICE_IN_KEY] = consts::DEFAULT_AUDIO_DEVICE_IN;

      if (config.find(consts::LEAD_CHANNEL_OUT_KEY) == config.end())
        config[consts::LEAD_CHANNEL_OUT_KEY] = consts::DEFAULT_LEAD_CHANNEL_OUT;

      if (config.find(consts::LEAD_CHANNEL_IN_KEY) == config.end())
        config[consts::LEAD_CHANNEL_IN_KEY] = consts::DEFAULT_LEAD_CHANNEL_IN;

    } else {
      config[consts::SOUND_OUTPUT_PATH_KEY] =
        al::File::conformPathToOS(userPath + consts::DEFAULT_SOUND_OUTPUT_PATH);

      config[consts::SAMPLE_RATE_KEY] = consts::SAMPLE_RATE;

      config[consts::BUFFER_SIZE_KEY] = consts::BUFFER_SIZE;

      // config[con::FONT_SCALE_KEY] = con::FONT_SCALE;

      config[consts::WINDOW_WIDTH_KEY] = consts::WINDOW_WIDTH;

      config[consts::WINDOW_HEIGHT_KEY] = consts::WINDOW_HEIGHT;

      config[consts::IS_FIRST_LAUNCH_KEY] = consts::IS_FIRST_LAUNCH;

      config[consts::DEFAULT_AUDIO_DEVICE_OUT_KEY] = consts::DEFAULT_AUDIO_DEVICE_OUT;

      config[consts::DEFAULT_AUDIO_DEVICE_IN_KEY] = consts::DEFAULT_AUDIO_DEVICE_IN;

      config[consts::LEAD_CHANNEL_OUT_KEY] = consts::DEFAULT_LEAD_CHANNEL_OUT;

      config[consts::LEAD_CHANNEL_IN_KEY] = consts::DEFAULT_LEAD_CHANNEL_IN;
    }
    std::ofstream file((userPath + configFile).c_str());
    if (file.is_open()) file << config;

    return false;
  }

  json jsonReadConfig() {
    json config;

    std::ifstream ifs(userPath + configFile);

    if (ifs.is_open()) config = json::parse(ifs);

    return config;
  }

  template <typename T>
  bool jsonWriteToConfig(T value, std::string key) {
    json config;

    std::ifstream ifs(userPath + configFile);

    if (ifs.is_open()) config = json::parse(ifs);

    config[key] = value;

    std::ofstream file((userPath + configFile).c_str());

    if (file.is_open()) {
      file << config;
      return true;
    } else {
      return false;
    }
  }

  json config;

  std::string getExecutablePath() {
#if _WIN32
    char *exePath;
    if (_get_pgmptr(&exePath) != 0) exePath = "";

#elif __linux__
    char exePath[PATH_MAX];
    ssize_t len = ::readlink("/proc/self/exe", exePath, sizeof(exePath));
    if (len == -1 || len == sizeof(exePath)) len = 0;
    exePath[len] = '\0';
#else
    char exePath[PATH_MAX];
    uint32_t len = sizeof(exePath);
    if (_NSGetExecutablePath(exePath, &len) != 0) {
      exePath[0] = '\0';  // buffer too small (!)
    } else {
      // resolve symlinks, ., .. if possible
      char *canonicalPath = realpath(exePath, NULL);
      if (canonicalPath != NULL) {
        strncpy(exePath, canonicalPath, len);
        free(canonicalPath);
      }
    }
#endif
    return std::string(exePath);
  }


  std::string getContentPath_OSX(std::string s) {
  char delim = '/';
  size_t counter = 0;
  size_t i = s.size() - 1;
  while (counter < 2) {
    if (s[i] == delim) counter++;
    i--;
  }
  return s.substr(0, i + 2);
}

  std::string getUserHomePath() {
  char homedir[PATH_MAX];
#ifdef _WIN32
  snprintf(homedir, sizeof(homedir), "%s%s", getenv("HOMEDRIVE"), getenv("HOMEPATH"));
#else
  snprintf(homedir, sizeof(homedir), "%s", getenv("HOME"));
#endif
  std::string result = strdup(homedir);
  return result;
}

  IMGUI_API bool InputText(const char *label, std::string *str, ImGuiInputTextFlags flags = 0,
                           ImGuiInputTextCallback callback = NULL, void *user_data = NULL);


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
  int currentbin, lastbin = 0;
};
