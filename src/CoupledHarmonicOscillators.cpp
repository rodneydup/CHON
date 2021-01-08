#include <array>

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

// NUMBER OF PARTICLES
#define nX (4)  // NUMBER OF PARTICLES on X Axis

struct MyApp : public App {
  MyApp() {
    for (int i = 0; i <= nX; i++) {
      // Create element
      auto *newSpring = new Spring;
      springs.push_back(newSpring);
      // Register its parameter bundle with the ControlGUI
      gui << newSpring->bundle;
    }
  }

  Mesh mesh;
  std::array<Particle, nX + 2> particle;  // create particles (plus 2 boundary particles)
  int picked = -1;                        // variable keeping track of which particle is selected
  std::array<double, nX + 1> k;           // Spring Constants
  double springLength = 50.0f / (nX + 1);
  bool freedom[3] = {1, 0, 0};
  Texture texBlur;  // blurring filter for graph
  int amCounter;    // AM incrementor

  Reverb<float> reverb;
  gam::NoiseWhite<> tick;
  bool drawGUI = 1;

  Parameter mAll{"Mass All", "physics", 1.0f, "", 1.0f, 100.0f};  // Master mass
  Parameter b{"Damping", "physics", 0.0f, "", 0.0f, 5.0f};        // damping
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

  void onInit() {  // Called on app start
    std::cout << "onInit()" << std::endl;
    reverb.bandwidth(0.9);  // Low-pass amount on input, in [0,1]
    reverb.damping(0.1);    // High-frequency damping, in [0,1]
    reverb.decay(0.15);     // Tail decay factor, in [0,1]
    reverb.diffusion(0.72, 0.69, 0.707, 0.71);
    client.open(port, addr);
    texBlur.filter(Texture::LINEAR);
    nav().pullBack(60);
    nav().pos(0, 7, 0);
    nav().setHome();
    xFree.registerChangeCallback([&](bool x) {
      if (!x) {
        for (int i = 1; i <= nX; i++) {
          particle[i].x(particle[i].equilibrium[0]);
          particle[i].velocity[0] = 0;
        }
      }
    });
    yFree.registerChangeCallback([&](bool y) {
      if (!y) {
        for (int i = 1; i <= nX; i++) {
          particle[i].y(particle[i].equilibrium[1]);
          particle[i].velocity[1] = 0;
        }
      }
    });
    zFree.registerChangeCallback([&](bool z) {
      if (!z) {
        for (int i = 1; i <= nX; i++) {
          particle[i].z(particle[i].equilibrium[2]);
          particle[i].velocity[2] = 0;
        }
      }
    });
    graphAxis.setElements({"x", "y", "z"});
    fmAxis.setElements({"x", "y", "z"});
    amAxis.setElements({"x", "y", "z"});
    bellAxis.setElements({"x", "y", "z"});
  }

  void onCreate() override {  // Called when graphics context is available
    std::cout << "onCreate()" << std::endl;

    imguiInit();

    addIcosphere(mesh, springLength / 5, 4);
    mesh.generateNormals();

    // const char *pentScale[25] = {"c3", "d3", "f3", "g3", "a3", "c4", "d4", "f4", "g4", "a4",
    // "c5", "d5", "f5", "g5", "a5", "c6", "d6", "f6", "g6", "a6", "c7", "d7", "f7", "g7",
    // "a7"};
    const char *majScale[35] = {"c3", "d3", "e3", "f3", "g3", "a3", "b3", "c4", "d4",
                                "e4", "f4", "g4", "a4", "b4", "c5", "d5", "e5", "f5",
                                "g5", "a5", "b5", "c6", "d6", "e6", "f6", "g6", "a6",
                                "b6", "c7", "d7", "e7", "f7", "g7", "a7", "b7"};
    int OTSeries[100];
    // float majScale[100];
    float majIntervals[7] = {1, 9.0 / 8, 5.0 / 4, 4.0 / 3, 3.0 / 2, 5.0 / 3, 15.0 / 8};

    for (int i = 0; i < 100; i++) {
      OTSeries[i] = 60 * (i + 1);
      // majScale[i] = 60 * pow(2, floor(((i + 1) % 50) / 7.0)) * majIntervals[i % 7];
    }

    for (int i = 0; i <= nX + 1; i++) {
      particle[i].particle.set(mesh);
      particle[i].equilibrium[0] = (i * springLength) - ((springLength * (nX + 1)) / 2);
      particle[i].equilibrium[1] = 0;
      particle[i].equilibrium[2] = 0;
      particle[i].particle.pose.setPos(particle[i].equilibrium);
      particle[i].graph.primitive(Mesh::LINE_STRIP);
      particle[i].oscillator.freq(100 * (i));
      particle[i].bell.freq(gam::scl::freq(majScale[i % (sizeof(majScale) / 8)]));
      // particle[i + 1].bell.freq(majScale[i % 100]);
      particle[i].mass = mAll;
      particle[i].amSmooth.setTime(40.0f);
      particle[i].fmSmooth.setTime(40.0f);
    }

    navControl().useMouse(false);
  }

  void onAnimate(double dt) override {  // Called once before drawing

    float w = width();
    float h = height();
    for (int i = 0; i <= nX; i++) {
      k[i] = springs[i]->k;
    }

    freedom[0] = xFree;
    freedom[1] = yFree;
    freedom[2] = zFree;

    if (!pause) {
      updateVelocities(particle, springLength, freedom, k, mAll, b, 60);

      for (int i = 1; i <= nX; i++) {  // animate stuff
        if (i == picked) {
          for (int j = 0; j < 3; j++) particle[i].velocity[j] = 0;
          if (!xFree) particle[i].x(particle[i].equilibrium[0]);
          if (!yFree) particle[i].y(particle[i].equilibrium[1]);
          if (!zFree) particle[i].z(particle[i].equilibrium[2]);
        } else {  // add velocities
          particle[i].addVelocity();
          // std::cout << "got here!" << std::endl;
        }

        particle[i].updateDisplacement();

        if (DrawGraph) {
          particle[i].graph.translate(-w / (60 * graphSpeed), 0,
                                      0);  // move previous graph data left
          if (particle[i].graph.vertices().size() > w)
            particle[i].graph.vertices().erase(
              particle[i].graph.vertices().begin());               // erase left hand graph data
          float graphY = h * particle[i].displacement[graphAxis];  // calculate new phase value
          if (graphAxis == 0)
            graphY /= springLength * 2;  // scale for longitudinal waves
          else
            graphY /= 20;  // scale for transverse waves
          graphY += (h * ((graphSpread * i / (nX + 1)) +
                          ((0.75 - (graphSpread / 2)))));  // add offset for each graph
          particle[i].graph.vertex(w, graphY, 0);          // update with new graph data on right
        }

        if (fm) {
          if (fmAxis == 0)
            particle[i].fmSmooth.setTarget(abs((particle[i].displacement[fmAxis] / springLength)));
          else
            particle[i].fmSmooth.setTarget(abs(particle[i].displacement[fmAxis] / 5));
        }

        if (am) {
          if (amAxis == 0)
            particle[i].amSmooth.setTarget(abs((particle[i].displacement[amAxis] / springLength)));
          else
            particle[i].amSmooth.setTarget(abs(particle[i].displacement[amAxis] / 5));
          if (particle[i].amSmooth.getTargetValue() > 1) particle[i].amSmooth.setTarget(1.0f);
        }

        for (int j = 0; j < 3; j++) {
          if (signbit(particle[i].prevDisplacement[j]) != signbit(particle[i].displacement[j]) &&
              abs(particle[i].prevDisplacement[j] - particle[i].displacement[j]) > 0.00001)
            particle[i].zeroTrigger[j] = 1;
        }
      }
    }
    for (int i = 1; i <= nX; i++) {
      // OSC
      if (oscOn) {
        particle[i].updateDisplacement();
        if (oscX)
          client.send("/displacementX", i, float(particle[i].displacement[0] / springLength));
        if (oscY) client.send("/displacementY", i, float(particle[i].displacement[1]));
        if (oscZ) client.send("/displacementZ", i, float(particle[i].displacement[2]));
      }
    }
  }

  void onDraw(Graphics &g) override {  // Draw function
    g.clear(0);
    // texBlur.resize(fbWidth(), fbHeight());
    // g.tint(0.9);
    // g.quadViewport(texBlur, -1, -1, 2, 2);
    // g.tint(1);
    g.depthTesting(true);
    g.lighting(true);
    g.color(0.5, 0.5, 0.5);
    g.polygonMode(Graphics::FILL);

    if (drawBoundaries) {
      g.color(1);
      particle[0].particle.drawMesh(g);
      particle[nX + 1].particle.drawMesh(g);
    }

    if (drawParticles) {
      for (int i = 1; i <= nX; i++) {
        g.color(HSV((float(i) / nX), 0.5, 1));
        if (!xFree) particle[i].x(particle[i].equilibrium[0]);
        if (!yFree) particle[i].y(particle[i].equilibrium[1]);
        if (!zFree) particle[i].z(particle[i].equilibrium[2]);
        particle[i].setPos(particle[i].x(), particle[i].y(), particle[i].z());
        particle[i].particle.drawMesh(g);
      }
      // texBlur.copyFrameBuffer();
    }

    if (DrawGraph) {
      g.pushCamera();
      g.camera(Viewpoint::ORTHO_FOR_2D);  // Ortho [0:width] x [0:height]
      g.lighting(false);
      for (int i = 1; i <= nX; i++) {
        g.color(HSV((float(i) / nX), 0.5, 1));
        g.draw(particle[i].graph);
      }
      g.popCamera();
    }

    // handle hidpi displays for imgui (mostly for linux and windows)
    // ImGui::GetIO().FontGlobalScale = getCurrentWindowScale();

    if (drawGUI) {
      imguiBeginFrame();

      ParameterGUI::beginPanel("Physics");
      gui.drawBundleGUI();
      ParameterGUI::drawParameter(&mAll);
      ParameterGUI::drawParameter(&b);
      ImGui::Text("Degrees of Freedom");
      ParameterGUI::drawParameterBool(&xFree);
      ImGui::SameLine();
      ParameterGUI::drawParameterBool(&yFree);
      ImGui::SameLine();
      ParameterGUI::drawParameterBool(&zFree);
      ParameterGUI::drawParameterBool(&pause);
      ParameterGUI::endPanel();

      ParameterGUI::beginPanel("Synthesis");
      ParameterGUI::drawParameterBool(&bellSynthOn);
      ImGui::SetNextItemWidth(35);
      ImGui::SameLine();
      ParameterGUI::drawMenu(&bellAxis);
      ParameterGUI::drawParameter(&bellVolume);
      ParameterGUI::drawParameterBool(&AdditiveSynthOn);
      ParameterGUI::drawParameter(&additiveVolume);
      ParameterGUI::drawParameter(&fundamental);
      ParameterGUI::drawParameterBool(&fm);
      ImGui::SameLine();
      ImGui::SetNextItemWidth(35);
      ParameterGUI::drawMenu(&fmAxis);
      ParameterGUI::drawParameter(&fmFreqMultiplier);
      ParameterGUI::drawParameter(&fmWidth);
      ParameterGUI::drawParameterBool(&am);
      ImGui::SameLine();
      ImGui::SetNextItemWidth(35);
      ParameterGUI::drawMenu(&amAxis);
      ParameterGUI::drawParameterBool(&reverbOn);
      ParameterGUI::endPanel();

      ParameterGUI::beginPanel("Display");
      ImGui::Text("Graph");
      ParameterGUI::drawParameterBool(&DrawGraph);
      ImGui::SameLine();
      ImGui::SetNextItemWidth(35);
      ParameterGUI::drawMenu(&graphAxis);
      ParameterGUI::drawParameter(&graphSpread);
      ParameterGUI::drawParameter(&graphSpeed);
      ImGui::Text("Particle System");
      ParameterGUI::drawParameterBool(&drawParticles);
      ParameterGUI::drawParameterBool(&drawBoundaries);
      ImGui::Text("Framerate %.3f", ImGui::GetIO().Framerate);
      ParameterGUI::endPanel();

      ParameterGUI::beginPanel("OSC");
      ParameterGUI::drawParameterBool(&oscOn);
      ParameterGUI::drawParameterBool(&oscX);
      ImGui::SameLine();
      ParameterGUI::drawParameterBool(&oscY);
      ImGui::SameLine();
      ParameterGUI::drawParameterBool(&oscZ);
      if (ImGui::InputText("IP Address", addr, IM_ARRAYSIZE(addr),
                           ImGuiInputTextFlags_EnterReturnsTrue))
        resetOSC();
      if (ImGui::InputInt("Port", &port, ImGuiInputTextFlags_EnterReturnsTrue)) resetOSC();
      ParameterGUI::endPanel();

      imguiEndFrame();

      imguiDraw();
    }
  }

  void onSound(AudioIOData &io) override {  // Audio callback
    while (io()) {
      double s = 0;
      // if (am) amCounter -= 1;

      if (AdditiveSynthOn) {
        for (int i = 1; i <= nX; i++) {  // add all oscillator samples to be sent to output
          if (i * fundamental < 20000) {
            if (fm) {
              particle[i].FM.freq((i + 3) * fundamental * fmFreqMultiplier);
              particle[i].fmSmooth.process();
              particle[i].oscillator.freq(
                (i * fundamental) + (particle[i].FM() * particle[i].fmSmooth.getCurrentValue() *
                                     fmWidth * 1000));  // set freq
            } else {
              particle[i].oscillator.freq((i + 3) * fundamental);  // set freq
            }
            double sampleToAdd = particle[i].oscillator() * additiveVolume;  // scale
            sampleToAdd = sampleToAdd / (double(i) + 1.0);                   // scale for harmonics
            if (am) {
              particle[i].amSmooth.process();  // progress smoothvalue
              sampleToAdd =
                sampleToAdd *
                particle[i].amSmooth.getCurrentValue();  // scale for amplitude modulation
            }
            s += sampleToAdd;
          }
        }
      }

      if (bellSynthOn) {
        for (int i = 1; i <= nX; i++) {
          if (particle[i].zeroTrigger[bellAxis]) {
            particle[i].bellEnv = 1;
          }
          if (particle[i].bellEnv > 0) particle[i].bellEnv -= 0.00005;
          float env = particle[i].getBellEnv();
          s += particle[i].bell() * 0.1 * bellVolume * env;
          particle[i].zeroTrigger[bellAxis] = 0;
        }
      }

      if (reverbOn) {
        // Compute two wet channels of reverberation
        float wet1, wet2;
        reverb(s, wet1, wet2);
        io.out(0) = wet1;
        io.out(1) = wet2;
      } else {
        io.out(0) = s;
        io.out(1) = s;
      }
    }
  }

  void onMessage(osc::Message &m) {  // OSC message callback
    m.print();
  }

  // helper functions to get scene ray from mouse events
  Vec3d unproject(Vec3d screenPos) {
    auto &g = graphics();
    auto mvp = g.projMatrix() * g.viewMatrix() * g.modelMatrix();
    Matrix4d invprojview = Matrix4d::inverse(mvp);
    Vec4d worldPos4 = invprojview.transform(screenPos);
    return worldPos4.sub<3>(0) / worldPos4.w;
  }

  Rayd getPickRay(int screenX, int screenY) {
    Rayd r;
    Vec3d screenPos;
    screenPos.x = (screenX * 1. / width()) * 2. - 1.;
    screenPos.y = ((height() - screenY) * 1. / height()) * 2. - 1.;
    screenPos.z = -1.;
    Vec3d worldPos = unproject(screenPos);
    r.origin().set(worldPos);

    screenPos.z = 1.;
    worldPos = unproject(screenPos);
    r.direction().set(worldPos);
    r.direction() -= r.origin();
    r.direction().normalize();
    return r;
  }

  bool onMouseMove(const Mouse &m) override {
    // make a ray from mouse location
    Rayd r = getPickRay(m.x(), m.y());
    for (int i = 1; i <= nX; i++) particle[i].particle.event(PickEvent(Point, r));
    return true;
  }
  bool onMouseDown(const Mouse &m) override {
    Rayd r = getPickRay(m.x(), m.y());
    for (int i = 1; i <= nX; i++) {
      particle[i].particle.event(PickEvent(Pick, r));
      if (particle[i].particle.selected == 1) picked = i;
    }
    return true;
  }
  bool onMouseDrag(const Mouse &m) override {
    Rayd r = getPickRay(m.x(), m.y());
    for (int i = 1; i <= nX; i++) particle[i].particle.event(PickEvent(Drag, r));
    return true;
  }
  bool onMouseUp(const Mouse &m) override {
    Rayd r = getPickRay(m.x(), m.y());
    for (int i = 1; i <= nX; i++) particle[i].particle.event(PickEvent(Unpick, r));
    picked = -1;
    return true;
  }

  bool onKeyDown(Keyboard const &k) override {
    switch (k.key()) {
      case 'p':
        pause = 1 - pause;
        break;
      case 'g':
        drawGUI = 1 - drawGUI;
        break;
      case 'r':
        nav().home();
        break;
      default:
        break;
    }
    return true;
  }

 private:
  std::vector<Spring *> springs;
  BundleGUIManager gui;
};

int main() {
  AudioDevice dev = AudioDevice::defaultOutput();
  dev.print();

  MyApp app;
  app.configureAudio(dev, dev.defaultSampleRate(), 1024, dev.channelsOutMax(), dev.channelsInMax());
  app.start();
  return 0;
}