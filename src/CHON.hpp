#include "./calculations.h"
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

  virtual void onAnimate(double dt) override;

  /**
   * @brief Draw rate processing of synth interface.
   */
  virtual void onDraw(al::Graphics &g) override;

  virtual Vec3d unproject(Vec3d screenPos);

  virtual Rayd getPickRay(int screenX, int screenY);

  virtual bool onMouseMove(const Mouse &m) override;

  virtual bool onMouseDown(const al::Mouse &m) override;

  virtual bool onMouseDrag(const Mouse &m) override;

  virtual bool onMouseUp(const Mouse &m) override;

  virtual bool onKeyDown(al::Keyboard const &k) override;

  virtual void chonReset();

  /** MIDI Stuff **/
  //   void initMIDI();
  //   void updateActiveMIDIParams(const al::MIDIMessage &m);

  /**
   * @brief Called everytime a MIDI message is sent.
   */
  //   virtual void onMIDIMessage(const al::MIDIMessage &m) override;

  //   virtual void onExit() override;

  Mesh mesh;           // mesh for drawing particles
  int nX = 4;          // number of particles on x axis
  int xParticles = 4;  // Parameter for x particles count
  int nY = 4;          // number of particles on y axis
  int yParticles = 4;  // Parameter for y particles count
  bool twoDimensions = false;
  std::vector<Spring *> springs;
  std::unique_ptr<BundleGUIManager> gui;

  std::vector<Particle> particle;  // create particles (plus 2 boundary particles)
  std::vector<double> k;           // Spring Constants
  double springLength;             // Spacing between particles
  bool freedom[3] = {1, 0, 0};     // Degrees of freedom for particle movement
  Texture texBlur;                 // blurring filter for graph
  int amCounter;                   // AM incrementor
  int picked = -1;                 // variable keeping track of which particle is selected
  std::mutex resetLock;

  Reverb<float> reverb;
  gam::NoiseWhite<> tick;
  bool drawGUI = 1;

  Parameter mAll{"Mass All", "physics", 1.0f, "", 1.0f, 100.0f};  // Master mass
  Parameter b{"Damping", "physics", 0.0f, "", 0.0f, 5.0f};        // damping
  ParameterBool xFree{"X axis", "Degrees of Freedom", 1};
  ParameterBool yFree{"Y axis", "Degrees of Freedom", 0};
  ParameterBool zFree{"Z axis", "Degrees of Freedom", 0};
  ParameterBool pause{"Pause (press p)", "physics", 0};           // Pause Simulation
  ParameterBool DrawGraph{"Draw Graph", "draw", 1};               // Toggle Drawing Graph
  Parameter graphSpread{"Spread", "draw", 0.0f, "", 0.0f, 1.0f};  // Graph spread
  Parameter graphSpeed{"Speed", "draw", 10.0f, "", 1.0f, 30.0f};  // Graph draw speed
  ParameterMenu graphAxis{"Graph Axis"};  // choose which axis displacement to draw
  ParameterBool drawParticles{"Draw Particles", "draw", 1};    // Toggle drawing particles
  ParameterBool drawBoundaries{"Draw Boundaries", "draw", 0};  // Toggle drawing boundaries

  ParameterBool AdditiveSynthOn{"Additive Synth On", "Synthesis", 0};  // Additive Synth toggle
  ParameterBool bellSynthOn{"Bell Synth On", "Synthesis", 0};          // Additive Synth toggle
  ParameterMenu bellAxis{"Bell Axis"};
  Parameter bellVolume{"Bell Volume", "Synthesis", 0.5f, "", 0.0f, 1.0f};  // Volume of bell synth
  Parameter additiveVolume{
    "Additive Volume", "Synthesis", 0.5f, "", 0.0f, 1.0f};  // Volume of bell synth
  ParameterBool reverbOn{"Reverb On", "Synthesis", 1};      // Reverb
  ParameterBool fm{"FM", "Synthesis", 0};                   // FM toggle
  Parameter fmFreqMultiplier{
    "FM Freq", "Synthesis", 1.5f,
    "",        0.0f,        2.0f};  // FM Frequency Multiplier (relative to carrier freq)
  ParameterMenu fmAxis{"FM Axis"};
  ParameterBool am{"AM", "Synthesis", 0};  // AM toggle
  ParameterMenu amAxis{"AM Axis"};
  Parameter fmWidth{"FM Width", "Synthesis", 2.0f, "", 0.1f, 5.0f};             // FM Width
  Parameter fundamental{"Fundamental", "Synthesis", 60.0f, "", 1.0f, 2000.0f};  // Fundamental

  const char *pentScale[25] = {"c3", "d3", "f3", "g3", "a3", "c4", "d4", "f4", "g4",
                               "a4", "c5", "d5", "f5", "g5", "a5", "c6", "d6", "f6",
                               "g6", "a6", "c7", "d7", "f7", "g7", "a7"};
  const char *majScale[35] = {"c3", "d3", "e3", "f3", "g3", "a3", "b3", "c4", "d4",
                              "e4", "f4", "g4", "a4", "b4", "c5", "d5", "e5", "f5",
                              "g5", "a5", "b5", "c6", "d6", "e6", "f6", "g6", "a6",
                              "b6", "c7", "d7", "e7", "f7", "g7", "a7", "b7"};
  int OTSeries[100];
  float majIntervals[7] = {1, 9.0 / 8, 5.0 / 4, 4.0 / 3, 3.0 / 2, 5.0 / 3, 15.0 / 8};

  int port = 16447;             // osc port
  char addr[10] = "127.0.0.1";  // ip address
  osc::Send client;             // create an osc client
  ParameterBool oscOn{"OSC On", "OSC", 0};
  ParameterBool oscX{"X disp", "OSC", 0};
  ParameterBool oscY{"Y disp", "OSC", 0};
  ParameterBool oscZ{"Z disp", "OSC", 0};

  void resetOSC() {
    client.open(port, addr);
    std::cout << "got here" << std::endl;
  }
};
