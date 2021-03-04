// To Do:
// - Add 2D back into the app
// - Allow more variables to be sent by OSC (velocity, potential energy)
// - add presets saving and loading
// - Work on GUI

#include "CHON.hpp"

#include <array>

void CHON::onInit() {  // Called on app start
  std::cout << "onInit()" << std::endl;

  gui = std::make_unique<BundleGUIManager>();

  for (int i = 0; i <= nX; i++) {
    // Create element
    auto *newSpring = new Spring;
    springs.push_back(newSpring);
    // Register its parameter bundle with the ControlGUI
    *gui << newSpring->bundle;
  }

  particle.resize(nX + 2);
  k.resize(nX + 1);
  springLength = 50.0f / (nX + 1);

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

void CHON::onCreate() {  // Called when graphics context is available
  std::cout << "onCreate()" << std::endl;

  imguiInit();

  addIcosphere(mesh, springLength / 5, 4);
  mesh.generateNormals();

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

void CHON::chonReset() {
  resetLock.lock();

  nX = xParticles;

  particle.clear();
  k.clear();
  particle.resize(nX + 2);
  k.resize(nX + 1);
  springLength = 50.0f / (nX + 1);

  mesh.reset();
  addIcosphere(mesh, springLength / 5, 4);
  mesh.generateNormals();

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
  gui.reset(new BundleGUIManager());
  springs.clear();
  for (int i = 0; i <= nX; i++) {
    // Create element
    auto *newSpring = new Spring;
    springs.push_back(newSpring);
    // Register its parameter bundle with the ControlGUI

    *gui << newSpring->bundle;
  }

  resetLock.unlock();
}

void CHON::onAnimate(double dt) {  // Called once before drawing

  if (nX != xParticles) chonReset();

  float w = width();
  float h = height();

  for (int i = 0; i < springs.size(); i++) {
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
      if (oscX) client.send("/displacementX", i, float(particle[i].displacement[0] / springLength));
      if (oscY) client.send("/displacementY", i, float(particle[i].displacement[1]));
      if (oscZ) client.send("/displacementZ", i, float(particle[i].displacement[2]));
    }
  }
}

void CHON::onDraw(Graphics &g) {  // Draw function
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
    ImGui::Text("Particle Count");
    ImGui::DragInt("x", &xParticles, 1.0f, 1, 100);
    // ImGui::SameLine();
    // ImGui::DragInt("y", &yParticles, 1.0f, 1, 100);
    ImGui::Checkbox("2D", &twoDimensions);
    gui->drawBundleGUI();
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

void CHON::onSound(AudioIOData &io) {  // Audio callback
  while (io()) {
    double s = 0;
    // if (am) amCounter -= 1;
    if (resetLock.try_lock()) {
      resetLock.unlock();
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
Vec3d CHON::unproject(Vec3d screenPos) {
  auto &g = graphics();
  auto mvp = g.projMatrix() * g.viewMatrix() * g.modelMatrix();
  Matrix4d invprojview = Matrix4d::inverse(mvp);
  Vec4d worldPos4 = invprojview.transform(screenPos);
  return worldPos4.sub<3>(0) / worldPos4.w;
}

Rayd CHON::getPickRay(int screenX, int screenY) {
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

bool CHON::onMouseMove(const Mouse &m) {
  // make a ray from mouse location
  Rayd r = getPickRay(m.x(), m.y());
  for (int i = 1; i <= nX; i++) particle[i].particle.event(PickEvent(Point, r));
  return true;
}
bool CHON::onMouseDown(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int i = 1; i <= nX; i++) {
    particle[i].particle.event(PickEvent(Pick, r));
    if (particle[i].particle.selected == 1) picked = i;
  }
  return true;
}
bool CHON::onMouseDrag(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int i = 1; i <= nX; i++) particle[i].particle.event(PickEvent(Drag, r));
  return true;
}
bool CHON::onMouseUp(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int i = 1; i <= nX; i++) particle[i].particle.event(PickEvent(Unpick, r));
  picked = -1;
  return true;
}

bool CHON::onKeyDown(Keyboard const &k) {
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

int main() {
  AudioDevice dev = AudioDevice::defaultOutput();
  dev.print();

  CHON app;
  app.configureAudio(dev, dev.defaultSampleRate(), 1024, dev.channelsOutMax(), dev.channelsInMax());
  app.start();
  return 0;
}