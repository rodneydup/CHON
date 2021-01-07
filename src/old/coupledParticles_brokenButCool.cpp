#include "al/app/al_App.hpp"
#include "Gamma/Oscillator.h"
#include "Gamma/Noise.h"
#include "al/graphics/al_Shapes.hpp"
#include "al/ui/al_ControlGUI.hpp"
#include "al/ui/al_Parameter.hpp"
#include "al/math/al_Ray.hpp"
#include "al/ui/al_BoundingBox.hpp"
#include "al/ui/al_Pickable.hpp"
#include "al/sound/al_Reverb.hpp"
#include "./calculations_brokenButCool.h"
#include <array>

using namespace al;

// NUMBER OF PARTICLES
#define n (200) // NUMBER OF PARTICLES
// NUMBER OF PARTICLES

struct MyApp : public App
{
  MyApp()
  {
    for (int i = 0; i <= n; i++)
    {
      // Create element
      auto *newSpring =
          new Spring; // This memory is not freed and it should be...
      springs.push_back(newSpring);
      // Register its parameter bundle with the ControlGUI
      gui << newSpring->bundle;
    }
  }

  Mesh mesh;
  std::array<Particle, n> particle;
  int picked = -1;            // variable keeping track of which particle is selected
  std::array<float, n + 1> k; // Spring Constants
  float springLength = 2.0f / (n + 1);
  bool freedom[3] = {1, 1, 0};
  Texture texBlur; // blurring filter for graph
  int amCounter;   //AM incrementor

  Reverb<float> reverb;
  gam::NoiseWhite<> tick;

  bool tickTrig = 0;

  Mesh boundaries;

  Parameter mAll{"Mass All", "physics", 1.0f, "", 0.1f, 10.0f}; // Master mass
  Parameter b{"Damping", "physics", 0.0f, "", 0.0f, 5.0f};      // damping
  ParameterBool pause{"Pause (press p)", "physics", 0};         // Pause Simulation

  ParameterBool DrawGraph{"Draw Graph", "draw", 0};           // Toggle Drawing Graph
  ParameterBool graphStack{"Stack Graph", "draw", 1};         // Graph lines stacked/adjacent toggle
  ParameterBool drawParticles{"Draw Particles", "draw", 1};   // Toggle drawing particles
  ParameterBool drawBoundaries{"Draw Boundaries", "draw", 0}; // Toggle drawing boundaries

  ParameterBool synthOn{"Additive Synth On", "Synthesis", 0};                    // Synth toggle
  Parameter volume{"Volume", "Synthesis", 0.5f, "", 0.0f, 1.0f};                 // Volume
  ParameterBool fm{"FM", "Synthesis", 0};                                        // FM toggle
  ParameterBool am{"AM", "Synthesis", 0};                                        // AM toggle
  Parameter fmWidth{"FM Width", "Synthesis", 2.0f, "", 0.1f, 5.0f};              // FM Width
  ParameterBool tickOn{"Tick", "Synthesis", 0};                                  // tick toggle
  Parameter fundamental{"Fundamental", "Synthesis", 100.0f, "", 20.0f, 2000.0f}; // Fundamental

  short port = 16447;             // osc port
  const char *addr = "127.0.0.1"; // ip address
  osc::Send client;               // create an osc client
  ParameterBool oscOn{"OSC On", "OSC", 0};
  ParameterString oscPort{"OSC Port", "OSC"};
  ParameterString oscIP{"OSC IP", "OSC"};

  void onInit()
  { // Called on app start
    std::cout << "onInit()" << std::endl;
    reverb.bandwidth(0.6); // Low-pass amount on input, in [0,1]
    reverb.damping(0.5);   // High-frequency damping, in [0,1]
    reverb.decay(0.1);     // Tail decay factor, in [0,1]
    reverb.diffusion(0.72, 0.69, 0.707, 0.71);
    client.open(port, addr);
    texBlur.filter(Texture::LINEAR);
    nav().pullBack(6);
  }

  void onCreate() override
  { // Called when graphics context is available
    std::cout << "onCreate()" << std::endl;

    imguiInit();

    addIcosphere(mesh, 0.2 / n, 4);
    mesh.generateNormals();

    for (int i = 0; i < n; i++)
    {
      particle[i].particle.set(mesh);
      particle[i].equilibrium[0] = ((2.0f * (i + 1) / (n + 1)) - 1);
      particle[i].equilibrium[1] = 0;
      particle[i].particle.pose.setPos(particle[i].equilibrium);
      particle[i].graph.primitive(Mesh::LINE_STRIP);
      particle[i].oscillator.freq(100 * (i + 1));
    }

    addCube(boundaries, 0, 0.2 / n);
    boundaries.translate(-2.4, 0, 0);
    addCube(boundaries, 0, 0.2 / n);
    boundaries.translate(1.2, 0, 0);
    boundaries.decompress();
    boundaries.generateNormals();

    navControl().useMouse(false);
  }

  void onAnimate(double dt) override
  { // Called once before drawing

    for (int i = 0; i <= n; i++)
    {
      k[i] = springs[i]->k;
    }

    if (!pause)
    {

      updateVelocity(particle, 1, springLength, freedom, k, mAll, b, 60);

      for (int i = 0; i < n; i++)
      {                                                                                        // animate stuff
        double yprev = particle[i].particle.pose.get().x() - ((2.0f * (i + 1) / (n + 1)) - 1); // get deviation from equilibrium before update

        if (i == picked)
        {
          particle[i].particle.pose.setPos(Vec3f(particle[i].particle.pose.get().x(), particle[i].particle.pose.get().y(), 0)); // Don't add velocity to ith particle if picked/dragged
        }
        else
        {
          particle[i].particle.pose.setPos(Vec3f(particle[i].particle.pose.get().x() + particle[i].velocity[0], particle[i].particle.pose.get().y() + particle[i].velocity[1], 0));
        }

        double y = particle[i].particle.pose.get().x() - ((2.0f * (i + 1) / (n + 1)) - 1); // get deviation after update
        y *= n / 2;                                                                        // scale according to number of particles

        if (DrawGraph)
        {
          particle[i].graph.translate(-0.01, 0, 0); // move previous graph data left
          if (particle[i].graph.vertices().size() > 220)
            particle[i].graph.vertices().erase(particle[i].graph.vertices().begin()); // erase left hand graph data
          float j = graphStack ? 0 : i;                                               // stack or spread graphs?
          particle[i].graph.vertex(1.2, ((j / n) * 2.0) + y + 0.5, 0);               // update with new graph data on right
        }

        fm ? particle[i].oscillator.freq(((i + 1) * fundamental) + (y * n * fmWidth * 100)) : particle[i].oscillator.freq((i + 1) * fundamental);

        if (am)
        {
          double amValPrev = particle[i].amVal;
          particle[i].amVal = abs(y);
          if (particle[i].amVal > 1)
            particle[i].amVal = 1;
          particle[i].amValIncrement = (particle[i].amVal - amValPrev) / (gam::sampleRate() / ImGui::GetIO().Framerate);
          amCounter = (gam::sampleRate() / ImGui::GetIO().Framerate);
        }

        if (tickOn && signbit(yprev) != signbit(y))
          tickTrig = 1;

        // OSC
        if (oscOn)
          client.send("/displacement", i, float(y));
      }
    }
  }

  void onDraw(Graphics &g) override
  { // Draw function
    g.clear(0, 0, 0);

    g.depthTesting(true);
    g.lighting(true);
    g.color(0.5, 0.5, 0.5);
    g.polygonMode(Graphics::FILL);

    if (drawBoundaries)
      g.draw(boundaries);

    if (drawParticles)
    {
      for (int i = 0; i < n; i++)
      {
        g.color(HSV((float(i) / n), 0.5, 1));
        particle[i].particle.pose.setPos(Vec3f(particle[i].particle.pose.get().pos()[0], particle[i].particle.pose.get().pos()[1], 0));
        particle[i].particle.drawMesh(g);
      }
    }
    if (DrawGraph)
    {
      g.lighting(false);
      for (int i = 0; i < n; i++)
      {
        g.color(HSV((float(i) / n), 0.5, 1));
        g.draw(particle[i].graph);
      }
    }
    texBlur.copyFrameBuffer();

    imguiBeginFrame();

    ParameterGUI::beginPanel("Physics");
    //  for (int i = 0; i <= n; i++) ParameterGUI::drawBundle(&springs[i].bundle);
    gui.drawBundleGUI();
    ParameterGUI::drawParameter(&mAll);
    ParameterGUI::drawParameter(&b);
    ParameterGUI::drawParameterBool(&pause);
    ParameterGUI::endPanel();

    ParameterGUI::beginPanel("Synthesis");
    // ParameterGUI::drawParameterBool(&tickOn);
    ParameterGUI::drawParameterBool(&synthOn);
    ParameterGUI::drawParameter(&volume);
    ParameterGUI::drawParameter(&fundamental);
    ParameterGUI::drawParameterBool(&fm);
    ParameterGUI::drawParameter(&fmWidth);
    ParameterGUI::drawParameterBool(&am);
    ParameterGUI::endPanel();

    ParameterGUI::beginPanel("Display");
    ImGui::Text("Graph");
    ParameterGUI::drawParameterBool(&DrawGraph);
    ParameterGUI::drawParameterBool(&graphStack);
    ImGui::Text("Particle System");
    ParameterGUI::drawParameterBool(&drawParticles);
    ParameterGUI::drawParameterBool(&drawBoundaries);
    ImGui::Text("Framerate %.3f",
                ImGui::GetIO().Framerate);
    ParameterGUI::endPanel();

    ParameterGUI::beginPanel("OSC");
    ParameterGUI::drawParameterBool(&oscOn);
    ParameterGUI::drawParameterString(&oscIP);
    ParameterGUI::drawParameterString(&oscPort);
    ParameterGUI::endPanel();

    imguiEndFrame();

    imguiDraw();
  }

  void onSound(AudioIOData &io) override
  { // Audio callback
    while (io())
    {
      float tickEnv;
      double s = 0;
      if (am)
        amCounter -= 1;

      if (synthOn)
      {
        for (int i = 0; i < n; i++)
        { // add all oscillator samples to be sent to output
          if (i * fundamental < 20000)
          {
            double sampleToAdd = particle[i].oscillator() * 0.1 * volume; // scale
            sampleToAdd = sampleToAdd / (float(i) + 1.0);                 // scale for harmonics
            if (am)
              sampleToAdd = sampleToAdd * (particle[i].amVal - (particle[i].amValIncrement * amCounter)); // scale for amplitude modulation
            s += sampleToAdd;
          }
        }
      }
      if (tickTrig)
      { // trigger envelope on noise to make a tick
        tickTrig = 0;
        tickEnv = 100;
      }
      s += (tickEnv / 200.0f) * tick(); // add tick to be sent to output
      if (tickEnv > 0)
        tickEnv -= 1; // decrement tick envelope

      if (tickEnv < 0)
        tickEnv = 0; // protect against negative values

      // Compute two wet channels of reverberation
      float wet1, wet2;
      reverb(s, wet1, wet2);

      io.out(0) = s;
      io.out(1) = s;
    }
  }

  void onMessage(osc::Message &m)
  { // OSC message callback
    m.print();
  }

  // helper functions to get scene ray from mouse events
  Vec3d unproject(Vec3d screenPos)
  {
    auto &g = graphics();
    auto mvp = g.projMatrix() * g.viewMatrix() * g.modelMatrix();
    Matrix4d invprojview = Matrix4d::inverse(mvp);
    Vec4d worldPos4 = invprojview.transform(screenPos);
    return worldPos4.sub<3>(0) / worldPos4.w;
  }

  Rayd getPickRay(int screenX, int screenY)
  {
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

  bool onMouseMove(const Mouse &m) override
  {
    // make a ray from mouse location
    Rayd r = getPickRay(m.x(), m.y());
    for (int i = 0; i < n; i++)
      particle[i].particle.event(PickEvent(Point, r));
    return true;
  }
  bool onMouseDown(const Mouse &m) override
  {
    Rayd r = getPickRay(m.x(), m.y());
    for (int i = 0; i < n; i++)
    {
      particle[i].particle.event(PickEvent(Pick, r));
      if (particle[i].particle.selected == 1)
        picked = i;
    }
    return true;
  }
  bool onMouseDrag(const Mouse &m) override
  {
    Rayd r = getPickRay(m.x(), m.y());
    for (int i = 0; i < n; i++)
      particle[i].particle.event(PickEvent(Drag, r));
    return true;
  }
  bool onMouseUp(const Mouse &m) override
  {
    Rayd r = getPickRay(m.x(), m.y());
    for (int i = 0; i < n; i++)
      particle[i].particle.event(PickEvent(Unpick, r));
    picked = -1;
    return true;
  }

  bool onKeyDown(Keyboard const &k) override
  {
    switch (k.key())
    {
    case 'p':
      pause = 1 - pause;
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

int main()
{
  MyApp app;
  app.configureAudio(44100, 512, 2, 0);

  app.start();
  return 0;
}