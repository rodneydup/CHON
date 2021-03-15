// To Do:
// - Allow more variables to be sent by OSC (velocity, potential energy)
// - add presets saving and loading
// - Work on GUI
// - Make the viewport smarter depending on window size and 2d vs 1d

#include "CHON.hpp"

#include <array>

void CHON::onInit() {  // Called on app start
  std::cout << "onInit()" << std::endl;

  xSpringGUI = std::make_unique<BundleGUIManager>();
  ySpringGUI = std::make_unique<BundleGUIManager>();

  for (int i = 0; i <= nX; i++) {
    // Create element
    auto *newSpring = new Spring{"X Springs"};
    xSprings.push_back(newSpring);
    // Register its parameter bundle with the ControlGUI
    *xSpringGUI << newSpring->bundle;
  }

  for (int i = 0; i <= nY; i++) {
    // Create element
    auto *newSpring = new Spring{"Y Springs"};
    ySprings.push_back(newSpring);
    // Register its parameter bundle with the ControlGUI
    *ySpringGUI << newSpring->bundle;
  }

  particle.resize(nX + 2);
  for (int y = 0; y <= nX + 1; y++) particle[y].resize(nY + 2);

  kX.resize(nX + 1);
  kY.resize(nY + 1);
  springLength = 1.0f / (nX + 1);

  reverb.bandwidth(0.9);  // Low-pass amount on input, in [0,1]
  reverb.damping(0.1);    // High-frequency damping, in [0,1]
  reverb.decay(0.15);     // Tail decay factor, in [0,1]
  reverb.diffusion(0.72, 0.69, 0.707, 0.71);
  client.open(port, addr);
  texBlur.filter(Texture::LINEAR);

  nav().pos(0, 0, 0);
  nav().pullBack(1.5);
  nav().setHome();

  xFree.registerChangeCallback([&](bool x) {
    if (!x) {
      for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
          particle[i][j].x(particle[i][j].equilibrium[0]);
          particle[i][j].velocity[0] = 0;
        }
      }
    }
  });
  yFree.registerChangeCallback([&](bool y) {
    if (!y) {
      for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
          particle[i][j].y(particle[i][j].equilibrium[1]);
          particle[i][j].velocity[1] = 0;
        }
      }
    }
  });
  zFree.registerChangeCallback([&](bool z) {
    if (!z) {
      for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
          particle[i][j].z(particle[i][j].equilibrium[2]);
          particle[i][j].velocity[2] = 0;
        }
      }
    }
  });
  graphAxis.setElements({"x", "y", "z"});
  fmAxis.setElements({"x", "y", "z"});
  amAxis.setElements({"x", "y", "z"});
  bellAxis.setElements({"x", "y", "z"});
  bellScale.setElements({"Pentatonic", "Major", "Chromatic", "Harmonics", "Bohlen-Pierce"});

  bellScale.registerChangeCallback([&](int tuning) {
    switch (tuning) {
      case 0:
        scale = &pentScale;
        break;
      case 1:
        scale = &majScale;
        break;
      case 2:
        scale = &chromScale;
        break;
      case 3:
        scale = &otSeries;
        break;
      case 4:
        scale = &bpScale;
        break;
      default:
        break;
    }
  });
}

void CHON::onCreate() {  // Called when graphics context is available
  std::cout << "onCreate()" << std::endl;

  imguiInit();

  bodyFont = ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF(&RobotoMedium_compressed_data,
                                                                  RobotoMedium_compressed_size, 16);
  titleFont = ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF(
    &RobotoMedium_compressed_data, RobotoMedium_compressed_size, 20);

  addIcosphere(mesh, springLength / 5, 4);
  mesh.generateNormals();
  for (int x = 0; x <= nX + 1; x++)
    for (int y = 0; y <= nY + 1; y++) {
      particle[x][y].particle.set(mesh);
      particle[x][y].equilibrium[0] = (x * springLength) - ((springLength * (nX + 1)) / 2);
      particle[x][y].equilibrium[1] = (y * springLength) - ((springLength * (nY + 1)) / 2);
      particle[x][y].equilibrium[2] = 0;
      particle[x][y].particle.pose.setPos(particle[x][y].equilibrium);
      particle[x][y].graph.primitive(Mesh::LINE_STRIP);
      particle[x][y].mass = mAll;
    }

  navControl().useMouse(false);
  // navControl().disable();
}

void CHON::chonReset() {
  std::cout << "Reset Particles" << std::endl;
  resetLock.lock();
  std::cout << yParticles << std::endl;

  // changing camera depending on if it's a 2D or 1D particle system
  if (nY == 1 && yParticles > 1) {
    nav().home();
    nav().pullBack(2 + pow(0.002 * (1900 - width()), 2));
    nav().turnF(-0.5);
    nav().faceToward(Vec3d{1, 1, -1});
  } else if (yParticles == 1) {
    nav().home();
    if (!drawGUI) {
      nav().pos(0, 0, 0);
      nav().pullBack(1.2 + pow(0.002 * (1900 - w), 2));
    }
    if (drawGUI) {
      nav().pos(-0.5 * (500 / w), 0, 0);
      nav().pullBack(1.5 + pow(0.002 * (1900 - w), 2));
    }
  } else if (nY > 1) {
    nav().pullBack(2 + pow(0.002 * (1900 - width()), 2));
    nav().faceToward(Vec3d{1, 1, -1});
  }

  // synchronize these variables
  nX = xParticles;
  nY = yParticles;

  // reset these arrays
  particle.clear();
  kX.clear();
  kY.clear();
  particle.resize(nX + 2);
  for (int y = 0; y <= nX + 1; y++) particle[y].resize(nY + 2);
  kX.resize(nX + 1);
  kY.resize(nY + 1);

  springLength = 1.0f / (std::max(nX, nY) + 1);

  mesh.reset();
  addIcosphere(mesh, springLength / 5, 4);
  mesh.generateNormals();

  for (int x = 0; x <= nX + 1; x++)
    for (int y = 0; y <= nY + 1; y++) {
      particle[x][y].particle.set(mesh);
      particle[x][y].equilibrium[0] = (x * springLength) - ((springLength * (nX + 1)) / 2);
      particle[x][y].equilibrium[1] = (y * springLength) - ((springLength * (nY + 1)) / 2);
      particle[x][y].equilibrium[2] = 0;
      particle[x][y].particle.pose.setPos(particle[x][y].equilibrium);
      particle[x][y].graph.primitive(Mesh::LINE_STRIP);
      particle[x][y].mass = mAll;
    }

  xSpringGUI.reset(new BundleGUIManager());
  ySpringGUI.reset(new BundleGUIManager());

  xSprings.clear();
  ySprings.clear();

  for (int i = 0; i <= nX; i++) {
    // Create element
    auto *newSpring = new Spring{"X Springs"};
    xSprings.push_back(newSpring);
    // Register its parameter bundle with the ControlGUI
    *xSpringGUI << newSpring->bundle;
  }
  for (int i = 0; i <= nY; i++) {
    // Create element
    auto *newSpring = new Spring{"Y Springs"};
    ySprings.push_back(newSpring);
    // Register its parameter bundle with the ControlGUI
    *ySpringGUI << newSpring->bundle;
  }

  resetLock.unlock();
}

void CHON::onAnimate(double dt) {  // Called once before drawing

  if (nX != xParticles) chonReset();

  if (nY != yParticles) chonReset();

  // Stuff to handle the camera when window size changes
  if (w != width()) {
    if (nY > 1) {
      if (!drawGUI) {
        nav().pos(0, 0, 0);
        nav().pullBack(2 + pow(0.002 * (1900 - width()), 2));
      }
      if (drawGUI) {
        nav().pos(-0.5 * (500 / width()), 0, 0);
        nav().pullBack(2.5 + pow(0.002 * (1900 - width()), 2));
      }
    } else {
      if (!drawGUI) {
        nav().pos(0, 0, 0);
        nav().pullBack(1.2 + pow(0.002 * (1900 - width()), 2));
      }
      if (drawGUI) {
        nav().pos(-0.5 * (500 / width()), 0, 0);
        nav().pullBack(1.5 + pow(0.002 * (1900 - width()), 2));
      }
    }
  }

  w = width();
  h = height();

  for (int i = 0; i < xSprings.size(); i++) {
    kX[i] = xSprings[i]->k;
  }

  for (int i = 0; i < ySprings.size(); i++) {
    kY[i] = ySprings[i]->k;
  }

  freedom[0] = xFree;
  freedom[1] = yFree;
  freedom[2] = zFree;

  if (!pause) {
    updateVelocities(particle, springLength, freedom, kX, kY, mAll, b, 60);

    for (int x = 1; x <= nX; x++)
      for (int y = 1; y <= nY; y++) {  // animate stuff
        if (((y - 1) * nX) + x == picked) {
          for (int j = 0; j < 3; j++) particle[x][y].velocity[j] = 0;
        } else {  // add velocities
          particle[x][y].addVelocity();
        }

        particle[x][y].updateDisplacement();

        if (DrawGraph) {
          particle[x][y].graph.translate(-w / (60 * graphSpeed), 0,
                                         0);  // move previous graph data left
          if (particle[x][y].graph.vertices().size() > w)
            particle[x][y].graph.vertices().erase(
              particle[x][y].graph.vertices().begin());               // erase left hand graph data
          float graphY = h * particle[x][y].displacement[graphAxis];  // calculate new phase value
          if (graphAxis == 0)
            graphY /= springLength * 2;  // scale for longitudinal waves
          else
            graphY /= 20;  // scale for transverse waves
          graphY += (h * ((graphSpread * x / (nX + 1)) +
                          ((0.75 - (graphSpread / 2)))));  // add offset for each graph
          particle[x][y].graph.vertex(w, graphY, 0);       // update with new graph data on right
        }

        if (fm) {
          if (nY > 1) {
            if (fmAxis == 0)
              particle[x][y].fmSmooth.setTarget(abs(particle[x][y].displacement[fmAxis]));
            else
              particle[x][y].fmSmooth.setTarget(abs(particle[x][y].displacement[fmAxis]));
          } else {
            if (fmAxis == 0)
              particle[x][y].fmSmooth.setTarget(
                abs(particle[x][y].displacement[fmAxis] / springLength));
            else
              particle[x][y].fmSmooth.setTarget(abs(particle[x][y].displacement[fmAxis] / 5));
          }
        }

        if (am) {
          if (amAxis == 0)
            particle[x][y].amSmooth.setTarget(
              abs((particle[x][y].displacement[amAxis] / springLength)));
          else
            particle[x][y].amSmooth.setTarget(abs(particle[x][y].displacement[amAxis] / 5));
          if (particle[x][y].amSmooth.getTargetValue() > 1) particle[x][y].amSmooth.setTarget(1.0f);
        }

        for (int j = 0; j < 3; j++) {
          if (signbit(particle[x][y].prevDisplacement[j]) !=
                signbit(particle[x][y].displacement[j]) &&
              abs(particle[x][y].prevDisplacement[j] - particle[x][y].displacement[j]) > 0.00001)
            particle[x][y].zeroTrigger[j] = 1;
        }

        // OSC
        if (oscOn) {
          particle[x][y].updateDisplacement();
          if (oscX)
            client.send("/displacementX", ((y - 1) * nX) + x,
                        float(particle[x][y].displacement[0] / springLength));
          if (oscY)
            client.send("/displacementY", ((y - 1) * nX) + x,
                        float(particle[x][y].displacement[1]));
          if (oscZ)
            client.send("/displacementZ", ((y - 1) * nX) + x,
                        float(particle[x][y].displacement[2]));
        }
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
    if (nY > 1)
      for (int x = 1; x <= nX; x++) {
        particle[x][0].particle.drawMesh(g);
        particle[x][nY + 1].particle.drawMesh(g);
      }
    for (int y = 1; y <= nY; y++) {
      particle[0][y].particle.drawMesh(g);
      particle[nX + 1][y].particle.drawMesh(g);
    }
  }

  if (drawParticles) {
    for (int x = 1; x <= nX; x++)
      for (int y = 1; y <= nY; y++) {
        g.color(HSV((float(x * y) / (nX * nY)), 0.5, 1));
        if (!xFree) particle[x][y].x(particle[x][y].equilibrium[0]);
        if (!yFree) particle[x][y].y(particle[x][y].equilibrium[1]);
        if (!zFree) particle[x][y].z(particle[x][y].equilibrium[2]);
        particle[x][y].setPos(particle[x][y].x(), particle[x][y].y(), particle[x][y].z());
        particle[x][y].particle.drawMesh(g);
      }
    // texBlur.copyFrameBuffer();
  }

  if (DrawGraph) {
    g.pushCamera();
    g.camera(Viewpoint::ORTHO_FOR_2D);  // Ortho [0:width] x [0:height]
    g.lighting(false);
    for (int x = 1; x <= nX; x++)
      for (int y = 1; y <= nY; y++) {
        g.color(HSV((float(x * y) / (nX * nY)), 0.5, 1));
        g.draw(particle[x][y].graph);
      }
    g.popCamera();
  }

  // handle hidpi displays for imgui (mostly for linux and windows)
  // ImGui::GetIO().FontGlobalScale = getCurrentWindowScale();

  if (drawGUI) {
    imguiBeginFrame();
    ImGui::PushFont(bodyFont);

    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Display", 0, 0, 350, 250, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
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
    ImGui::Text("%.2f", nav().z());
    ImGui::PopFont();
    ParameterGUI::endPanel();

    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Physics", 0, 250, 350, 300, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    ImGui::Text("Particle Count");
    ImGui::PushItemWidth(100);
    ImGui::InputInt("x", &xParticles);
    if (xParticles * yParticles > 100)
      xParticles = int(100 / yParticles);
    else if (xParticles < 1)
      xParticles = 1;
    ImGui::SameLine();
    ImGui::InputInt("y", &yParticles);
    if (xParticles * yParticles > 100)
      yParticles = int(100 / xParticles);
    else if (yParticles < 1)
      yParticles = 1;
    ImGui::PopItemWidth();
    ImGui::PushItemWidth(200);
    xSpringGUI->drawBundleGUI();
    if (nY > 1) ySpringGUI->drawBundleGUI();
    ImGui::PopItemWidth();
    ParameterGUI::drawParameter(&mAll);
    ParameterGUI::drawParameter(&b);
    ImGui::Text("Degrees of Freedom");
    ParameterGUI::drawParameterBool(&xFree);
    ImGui::SameLine();
    ParameterGUI::drawParameterBool(&yFree);
    ImGui::SameLine();
    ParameterGUI::drawParameterBool(&zFree);
    ParameterGUI::drawParameterBool(&pause);
    ImGui::PopFont();
    ParameterGUI::endPanel();

    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Synthesis", 0, 550, 350, 275, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    ParameterGUI::drawParameterBool(&bellSynthOn);
    if (bellSynthOn) {
      ImGui::SetNextItemWidth(35);
      ImGui::SameLine();
      ParameterGUI::drawMenu(&bellAxis);
      ParameterGUI::drawMenu(&bellScale);
      ParameterGUI::drawParameter(&bellRoot);
      ParameterGUI::drawParameter(&bellVolume);
    }
    ParameterGUI::drawParameterBool(&AdditiveSynthOn);
    if (AdditiveSynthOn) {
      ParameterGUI::drawParameter(&additiveVolume);
      ParameterGUI::drawParameter(&additiveRoot);
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
    }
    ParameterGUI::drawParameterBool(&reverbOn);
    ImGui::PopFont();
    ParameterGUI::endPanel();

    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("OSC", 0, 825, 350, 150, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
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
    ImGui::PopFont();
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
        for (int x = 1; x <= nX; x++)
          for (int y = 1; y <= nY; y++) {  // add all oscillator samples to be sent to output
            if (x * additiveRoot < 20000) {
              if (fm) {
                particle[x][y].FM.freq(x * y * additiveRoot * fmFreqMultiplier);
                particle[x][y].fmSmooth.process();
                particle[x][y].oscillator.freq(
                  (x * y * additiveRoot) +
                  (particle[x][y].FM() * particle[x][y].fmSmooth.getCurrentValue() * fmWidth *
                   1000));  // set freq
              } else {
                particle[x][y].oscillator.freq(x * y * additiveRoot);  // set freq
              }
              double sampleToAdd = particle[x][y].oscillator() * additiveVolume;  // scale
              sampleToAdd = sampleToAdd / ((x * y) + 1.0);  // scale for harmonics
              if (am) {
                particle[x][y].amSmooth.process();  // progress smoothvalue
                sampleToAdd =
                  sampleToAdd *
                  particle[x][y].amSmooth.getCurrentValue();  // scale for amplitude modulation
              }
              s += sampleToAdd;
            }
          }
      }

      if (bellSynthOn) {
        for (int x = 1; x <= nX; x++)
          for (int y = 1; y <= nY; y++) {
            if (particle[x][y].zeroTrigger[bellAxis]) {
              if ((scale->at((x * y) % (scale->size())) * bellRoot) < 20000)
                particle[x][y].bell.freq(scale->at((x * y) % (scale->size())) * bellRoot);
              else
                particle[x][y].bell.freq(0);
              particle[x][y].bellEnv = 1;
            }
            if (particle[x][y].bellEnv > 0) particle[x][y].bellEnv -= 0.00005;
            s += particle[x][y].bell() * 0.2 * bellVolume * particle[x][y].getBellEnv();
            particle[x][y].zeroTrigger[bellAxis] = 0;
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
  for (int x = 1; x <= nX; x++)
    for (int y = 1; y <= nY; y++) particle[x][y].particle.event(PickEvent(Point, r));
  return true;
}
bool CHON::onMouseDown(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int x = 1; x <= nX; x++)
    for (int y = 1; y <= nY; y++) {
      particle[x][y].particle.event(PickEvent(Pick, r));
      if (particle[x][y].particle.selected == 1) picked = ((y - 1) * nX) + x;
    }
  return true;
}
bool CHON::onMouseDrag(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int x = 1; x <= nX; x++)
    for (int y = 1; y <= nY; y++) particle[x][y].particle.event(PickEvent(Drag, r));
  return true;
}
bool CHON::onMouseUp(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int x = 1; x <= nX; x++)
    for (int y = 1; y <= nY; y++) particle[x][y].particle.event(PickEvent(Unpick, r));
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
      if (nY > 1) {
        if (!drawGUI) {
          nav().pos(0, 0, 0);
          nav().pullBack(1.2 + pow(0.002 * (1900 - w), 2));
        }
        if (drawGUI) {
          nav().pos(-0.5 * (500 / w), 0, 0);
          nav().pullBack(1.5 + pow(0.002 * (1900 - w), 2));
        }
      } else {
        if (!drawGUI) {
          nav().pos(0, 0, 0);
          nav().pullBack(1.2 + pow(0.002 * (1900 - w), 2));
        }
        if (drawGUI) {
          nav().pos(-0.5 * (500 / w), 0, 0);
          nav().pullBack(1.5 + pow(0.002 * (1900 - w), 2));
        }
      }
      break;
    case 'r':
      // nav().home();
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