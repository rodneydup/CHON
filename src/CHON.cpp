// To Do:
// FEATURES:
//    - Allow more variables to be sent by OSC (velocity, potential energy)
//    - add presets saving and loading
//    - Implement "kuramoto mode"
//    - Implement collisions
//    - Accept OSC in to control some parameters

#include "CHON.hpp"

#include <array>

void CHON::onInit() {  // Called on app start
  std::cout << "onInit()" << std::endl;
  title("CHON");

  reverb.bandwidth(0.9);  // Low-pass amount on input, in [0,1]
  reverb.damping(0.1);    // High-frequency damping, in [0,1]
  reverb.decay(0.15);     // Tail decay factor, in [0,1]
  reverb.diffusion(0.72, 0.69, 0.707, 0.71);
  client.open(port, addr);
  texBlur.filter(Texture::LINEAR);

  nav().pos(0, 0.2, 0);
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
  reverbTail.registerChangeCallback([&](float tail) { reverb.decay(tail); });

  graphAxis.setElements({"x", "y", "z"});
  fmAxis.setElements({"x", "y", "z"});
  amAxis.setElements({"x", "y", "z"});
  bellAxis.setElements({"x", "y", "z"});
  bellScale.setElements({"Pentatonic", "Major", "Chromatic", "Harmonics", "Bohlen-Pierce"});
  driveAxisLeft.setElements({"x", "y", "z"});
  driveAxisRight.setElements({"x", "y", "z"});
  inputMode.setElements({"Peak", "RMS", "Frequency"});

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
  audioIO().setStreamName("CHON");
}

void CHON::onCreate() {  // Called when graphics context is available
  std::cout << "onCreate()" << std::endl;

  imguiInit();
  ImGuiIO &io = ImGui::GetIO();
  io.IniFilename = NULL;

  bodyFont = ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF(&RobotoMedium_compressed_data,
                                                                  RobotoMedium_compressed_size, 16);
  titleFont = ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF(
    &RobotoMedium_compressed_data, RobotoMedium_compressed_size, 20);

  chonReset();

  navControl().useMouse(false);
  navControl().disable();
}

void CHON::chonReset() {
  std::cout << "Reset Particle network to " << xParticles << " by " << yParticles << std::endl;
  resetLock.lock();

  // changing camera depending on if it's a 2D or 1D particle system
  if (nY == 1 && yParticles > 1) {
    onResize(width(), height());
  } else if (yParticles == 1) {
    onResize(width(), height());
  }

  // synchronize these variables
  nX = xParticles;
  nY = yParticles;

  // reset these arrays
  particle.clear();
  particle.resize(nX + 2);
  for (int y = 0; y <= nX + 1; y++) particle[y].resize(nY + 2);

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
      particle[x][y].oscName[0] = "/dispX/" + std::to_string(((y - 1) * nX) + x);
      particle[x][y].oscName[1] = "/dispY/" + std::to_string(((y - 1) * nX) + x);
      particle[x][y].oscName[2] = "/dispZ/" + std::to_string(((y - 1) * nX) + x);
      particle[x][y].oscName[3] = "/pan/" + std::to_string(((y - 1) * nX) + x);
    }

  xSpringGUI.reset(new ChonBundle(1));
  ySpringGUI.reset(new ChonBundle(1));

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

  driveParticleXLeft.max(nX);
  driveParticleYLeft.max(nY);
  driveParticleXRight.max(nX);
  driveParticleYRight.max(nY);

  fftDivision = log(fftBuffer.size()) / (nX * nY);
  client.send("/Nparticles/", nX * nY);
  resetLock.unlock();
}

void CHON::onAnimate(double dt) {  // Called once before drawing

  if (nX != xParticles) chonReset();

  if (nY != yParticles) chonReset();

  w = width();
  h = height();

  freedom[0] = xFree;
  freedom[1] = yFree;
  freedom[2] = zFree;

  if (!pause) {
    if (inputOn) {
      switch (inputMode) {
        case 0:  // Peak
          driveForceLeft = inLeft.at(inLeft.getTail());
          if (driveForceLeft > inputThreshold)
            particle[driveParticleXLeft][driveParticleYLeft].acceleration[driveAxisLeft] +=
              driveForceLeft * inputScale;
          if (driveStereoSplit) {
            driveForceRight = inRight.at(inRight.getTail());
            if (driveForceRight > inputThreshold)
              particle[driveParticleXRight][driveParticleYRight].acceleration[driveAxisRight] +=
                driveForceRight * inputScale;
          }
          break;
        case 1:  // RMS
          driveForceLeft = inLeft.getRMS(rmsSize);
          if (driveForceLeft > inputThreshold)
            particle[driveParticleXLeft][driveParticleYLeft].acceleration[driveAxisLeft] +=
              driveForceLeft * inputScale;
          if (driveStereoSplit) {
            driveForceRight = inRight.getRMS(rmsSize);
            if (driveForceRight > inputThreshold)
              particle[driveParticleXRight][driveParticleYRight].acceleration[driveAxisRight] +=
                driveForceRight * inputScale;
          }
          break;
        case 2:  // Frequency
          particle[1][1].acceleration[driveAxisLeft] += fftBuffer[0] * inputScale;
          for (int i = 1; i < fftBuffer.size(); i++) {
            fftIterator = floor((log(i) / fftDivision) + 1);
            particle[fftIterator % (nX + 1)][ceil(fftIterator / float(nX + 1))]
              .acceleration[driveAxisLeft] += fftBuffer[i] * inputScale;
          }
          fftIterator = 0;
          break;
        default:
          break;
      }
    }

    updateVelocities(particle, springLength, freedom, xSprings, ySprings, mAll, b, dt);
    for (int x = 1; x <= nX; x++)
      for (int y = 1; y <= nY; y++) {  // animate stuff
        if (((y - 1) * nX) + x == picked) {
          for (int j = 0; j < 3; j++) particle[x][y].velocity[j] = 0;
        } else {  // increment positions according to velocity
          particle[x][y].addVelocity();
        }

        particle[x][y].updateDisplacement();

        if (DrawGraph) {  // updating the 2d displacement graph
          particle[x][y].graph.translate(-w / (60 * graphSpeed), 0,
                                         0);  // move previous graph data left
          if (particle[x][y].graph.vertices().size() > w)
            particle[x][y].graph.vertices().erase(
              particle[x][y].graph.vertices().begin());               // erase left-edge graph data
          float graphY = h * particle[x][y].displacement[graphAxis];  // calculate new phase value
          graphY /= springLength * 2;
          graphY += (h * ((graphSpread * x / (nX + 1)) +
                          ((0.75 - (graphSpread / 2)))));  // add offset for each graph
          particle[x][y].graph.vertex(w, graphY, 0);       // update with new graph data on right
        }

        if (fm) {  // copy displacement values to FM synthesis engine
          particle[x][y].fmSmooth.setTarget(particle[x][y].displacement[fmAxis] / springLength);
        }

        if (am) {  // copy displacement values to AM synthesis engine
          float val = particle[x][y].displacement[amAxis] / springLength;
          val > 1 ? 1 : val;
          val < -1 ? -1 : val;
          particle[x][y].amSmooth.setTarget(val);
        }

        if (stereoOn) {  // set the particle's pan position
          particle[x][y].panSmooth.setTarget(particle[x][y].x() + 0.5);
        }

        for (int j = 0; j < 3; j++) {  // check if zero crossing
          if (signbit(particle[x][y].prevDisplacement[j]) !=
                signbit(particle[x][y].displacement[j]) &&
              abs(particle[x][y].prevDisplacement[j] - particle[x][y].displacement[j]) > 0.00001) {
            particle[x][y].zeroTrigger[j] = 1;
          }
        }

        // OSC
        if (oscOn) {
          if (oscX)
            client.send(particle[x][y].oscName[0],
                        float(particle[x][y].displacement[0] / springLength));
          if (oscY) client.send(particle[x][y].oscName[1], float(particle[x][y].displacement[1]));
          if (oscZ) client.send(particle[x][y].oscName[2], float(particle[x][y].displacement[2]));
          if (oscPan) client.send(particle[x][y].oscName[3], float(particle[x][y].x()));
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
  g.polygonFill();

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

    int yposition = 0;
    ImGui::SetNextWindowCollapsed(1, ImGuiCond_Once);
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Display", 0, yposition, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    ParameterGUI::drawParameterBool(&DrawGraph);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(35);
    ParameterGUI::drawMenu(&graphAxis);
    ParameterGUI::drawParameter(&graphSpread);
    ParameterGUI::drawParameter(&graphSpeed);
    ParameterGUI::drawParameterBool(&drawParticles);
    ImGui::SameLine();
    ParameterGUI::drawParameterBool(&drawBoundaries);
    ImGui::Text("Framerate %.3f", ImGui::GetIO().Framerate);
    ImGui::PopFont();
    yposition += ImGui::GetWindowHeight();
    ParameterGUI::endPanel();

    ImGui::SetNextWindowCollapsed(1, ImGuiCond_Once);
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Physics", 0, yposition, flags);
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
    ImGui::Separator();
    ImGui::PushItemWidth(200);
    xSpringGUI->drawBundleGUI();
    if (nY > 1) ySpringGUI->drawBundleGUI();
    ImGui::PopItemWidth();
    ParameterGUI::drawParameter(&mAll);
    ParameterGUI::drawParameter(&b);
    ImGui::Separator();
    ImGui::Text("Degrees of Freedom");
    ParameterGUI::drawParameterBool(&xFree);
    ImGui::SameLine();
    ParameterGUI::drawParameterBool(&yFree);
    ImGui::SameLine();
    ParameterGUI::drawParameterBool(&zFree);
    ParameterGUI::drawParameterBool(&pause);
    ImGui::PopFont();
    yposition += ImGui::GetWindowHeight();
    ParameterGUI::endPanel();

    ImGui::SetNextWindowCollapsed(1, ImGuiCond_Once);
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Synthesis", 0, yposition, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    ParameterGUI::drawParameterBool(&stereoOn);
    if (ImGui::CollapsingHeader("Bell Synth")) {
      ParameterGUI::drawParameterBool(&bellSynthOn);
      ImGui::SetNextItemWidth(35);
      ImGui::SameLine();
      ParameterGUI::drawMenu(&bellAxis);
      ParameterGUI::drawMenu(&bellScale);
      ParameterGUI::drawParameter(&bellRoot);
      ParameterGUI::drawParameter(&bellVolume);
    }
    if (ImGui::CollapsingHeader("Additive Synth")) {
      ParameterGUI::drawParameterBool(&additiveSynthOn);
      ParameterGUI::drawParameter(&additiveRoot);
      ParameterGUI::drawParameter(&additiveVolume);
      if (ImGui::CollapsingHeader("FM")) {
        ParameterGUI::drawParameterBool(&fm);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(35);
        ParameterGUI::drawMenu(&fmAxis);
        ParameterGUI::drawParameter(&fmFreq);
        ParameterGUI::drawParameter(&fmWidth);
      }
      if (ImGui::CollapsingHeader("AM")) {
        ParameterGUI::drawParameterBool(&am);
        {
          ImGui::SameLine();
          ImGui::SetNextItemWidth(35);
          ParameterGUI::drawMenu(&amAxis);
        }
      }
    }
    if (ImGui::CollapsingHeader("Reverb")) {
      ParameterGUI::drawParameterBool(&reverbOn);
      ParameterGUI::drawParameter(&reverbTail);
    }
    ImGui::PopFont();
    yposition += ImGui::GetWindowHeight();
    ParameterGUI::endPanel();

    ImGui::SetNextWindowCollapsed(1, ImGuiCond_Once);
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("Audio", 0, yposition, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    if (ImGui::CollapsingHeader("Audio IO Settings")) {
      drawAudioIO(&audioIO());
    }
    if (ImGui::CollapsingHeader("Input")) {
      ParameterGUI::drawParameterBool(&inputOn);
      ParameterGUI::drawMenu(&inputMode);
      if (inputMode != 2) {
        ParameterGUI::drawParameterBool(&driveStereoSplit);
        if (driveStereoSplit) ImGui::Text("Left Channel");
        ParameterGUI::drawParameterInt(&driveParticleXLeft, "");
        ParameterGUI::drawParameterInt(&driveParticleYLeft, "");
      }
      ParameterGUI::drawMenu(&driveAxisLeft);
      if (driveStereoSplit && inputMode != 2) {
        ImGui::Text("Right Channel");
        ParameterGUI::drawParameterInt(&driveParticleXRight, "");
        ParameterGUI::drawParameterInt(&driveParticleYRight, "");
        ParameterGUI::drawMenu(&driveAxisRight);
      }
      ParameterGUI::drawParameter(&inputScale);
      if (inputMode != 2) ParameterGUI::drawParameter(&inputThreshold);
      if (inputMode == 1) ParameterGUI::drawParameterInt(&rmsSize, "");
    }
    yposition += ImGui::GetWindowHeight();
    ImGui::PopFont();
    ParameterGUI::endPanel();

    ImGui::SetNextWindowCollapsed(1, ImGuiCond_Once);
    ImGui::PushFont(titleFont);
    ParameterGUI::beginPanel("OSC", 0, yposition, flags);
    ImGui::PopFont();
    ImGui::PushFont(bodyFont);
    ParameterGUI::drawParameterBool(&oscOn);
    ImGui::Text("Sending:");
    ParameterGUI::drawParameterBool(&oscX);
    ImGui::SameLine();
    ParameterGUI::drawParameterBool(&oscY);
    ImGui::SameLine();
    ParameterGUI::drawParameterBool(&oscZ);
    ParameterGUI::drawParameterBool(&oscPan);
    if (ImGui::InputText("IP Address", addr, IM_ARRAYSIZE(addr),
                         ImGuiInputTextFlags_EnterReturnsTrue))
      resetOSC();
    if (ImGui::InputInt("Port", &port, ImGuiInputTextFlags_EnterReturnsTrue)) resetOSC();
    if (ImGui::CollapsingHeader("Message Syntax")) {
      ImGui::Text(
        "OSC messages are formatted as follows:\n"
        "FirstArg/SecondArg Value\n\n"
        "FirstArg:\n"
        "dispX, dispY, or dispZ \n\n"
        "SecondArg:\n"
        "Integer representing the particle number");
    }
    ImGui::PopFont();

    ParameterGUI::endPanel();

    imguiEndFrame();

    imguiDraw();
  }
}

void CHON::onSound(AudioIOData &io) {  // Audio callback
  while (io()) {
    double s[2] = {0, 0};

    if (resetLock.try_lock()) {
      if (stereoOn) {
        for (int x = 1; x <= nX; x++)
          for (int y = 1; y <= nY; y++) {
            particle[x][y].panSmooth.process();
          }
      }

      if (additiveSynthOn) {
        for (int x = 1; x <= nX; x++)
          for (int y = 1; y <= nY; y++) {  // add all oscillator samples to be sent to output
            if (x * additiveRoot < 20000) {
              if (fm) {
                particle[x][y].FM.freq(x * y * additiveRoot * fmFreq);
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
              if (stereoOn) {
                s[0] += (1 - particle[x][y].panSmooth.getCurrentValue()) * sampleToAdd;
                s[1] += particle[x][y].panSmooth.getCurrentValue() * sampleToAdd;
              } else {
                s[0] += sampleToAdd;
                s[1] += sampleToAdd;
              }
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
            double sampleToAdd =
              particle[x][y].bell() * 0.2 * bellVolume * particle[x][y].getBellEnv();
            if (stereoOn) {
              s[0] += (1 - particle[x][y].panSmooth.getCurrentValue()) * sampleToAdd;
              s[1] += particle[x][y].panSmooth.getCurrentValue() * sampleToAdd;
            } else {
              s[0] += sampleToAdd;
              s[1] += sampleToAdd;
            }

            particle[x][y].zeroTrigger[bellAxis] = 0;
          }
      }
      resetLock.unlock();
    }

    if (reverbOn) {
      // Compute two wet channels of reverberation
      float wet1, wet2, throwAway;
      reverb(s[0], wet1, throwAway);
      reverb(s[1], throwAway, wet2);
      io.out(0) = wet1;
      io.out(1) = wet2;
    } else {
      io.out(0) = s[0];
      io.out(1) = s[1];
    }

    // Handle audio input
    if (inputOn) {
      switch (inputMode) {
        case 0:  // Peak
          if (driveStereoSplit) {
            inLeft.push_back(io.in(0));
            inRight.push_back(io.in(1));
          } else {
            inLeft.push_back((io.in(0) + io.in(1)) / 2);
          }
          break;
        case 1:  // RMS
          if (driveStereoSplit) {
            inLeft.push_back(io.in(0));
            inRight.push_back(io.in(1));
          } else {
            inLeft.push_back((io.in(0) + io.in(1)) / 2);
          }
          break;
        case 2:  // Frequency
          if (stft((io.in(0) + io.in(1)) / 2)) {
            for (int i = 0; i < stft.numBins(); i++) {
              fftBuffer[i] = abs(stft.bin(i)[0]);
            }
          }
          break;

        default:
          break;
      }
    }
  }
}

void CHON::onResize(int width, int height) {
  if (!drawGUI) {
    nav().home();
    if (yParticles > 1) {
      nav().pullBack(2.2 + pow(0.002 * (1900 - width), 2));
      nav().view(0, 0, -0.5);
      nav().faceToward(Vec3d{1, 1, -1});
    } else
      nav().pullBack(1.2 + pow(0.002 * (1900 - width), 2));
  } else {
    nav().home();
    nav().pos(-0.5 * (500 / float(width)), 0.2, 0);
    if (yParticles > 1) {
      nav().pullBack(2.7 + pow(0.002 * (1900 - width), 2));
      nav().view(0, 0, -0.5);
      nav().faceToward(Vec3d{1, 1, -1});
    } else
      nav().pullBack(1.5 + pow(0.002 * (1900 - width), 2));
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
      return false;
    case 'v':
      onResize(width(), height());
      return false;
    case 'g':
      drawGUI = 1 - drawGUI;
      onResize(width(), height());
      return false;
    case 'r':
      srand(std::time(0));
      for (int x = 1; x <= nX; x++)
        for (int y = 1; y <= nY; y++) {
          particle[x][y].velocity[0] += (float(rand() - (RAND_MAX / 2)) / RAND_MAX) / 5;
          particle[x][y].velocity[1] += (float(rand() - (RAND_MAX / 2)) / RAND_MAX) / 5;
          particle[x][y].velocity[2] += (float(rand() - (RAND_MAX / 2)) / RAND_MAX) / 5;
        }
      return false;
    case Keyboard::UP:
      nav().moveU(0.02);
      return false;
    case Keyboard::DOWN:
      nav().moveU(-0.02);
      return false;
    case Keyboard::LEFT:
      nav().spinU(0.02);
      return false;
    case Keyboard::RIGHT:
      nav().spinU(-0.02);
      return false;
    default:
      break;
  }
  return true;
}

bool CHON::onKeyUp(Keyboard const &k) {
  switch (k.key()) {
    case Keyboard::UP:
      nav().moveU(0);
      return false;
    case Keyboard::DOWN:
      nav().moveU(0);
      return false;
    case Keyboard::LEFT:
      nav().spinU(0);
      return false;
    case Keyboard::RIGHT:
      nav().spinU(0);
    default:
      break;
  }
  return true;
}

void CHON::drawAudioIO(AudioIO *io) {
  struct AudioIOState {
    int currentSr = 1;
    int currentBufSize = 3;
    int currentDeviceOut = 0;
    int currentDeviceIn = 0;
    int currentOut = 1;
    int currentIn = 1;
    int currentMaxOut;
    int currentMaxIn;
    std::vector<std::string> devices;
  };

  auto updateOutDevices = [&](AudioIOState &state) {
    state.devices.clear();
    int numDevices = AudioDevice::numDevices();
    int dev_out_index = 0;
    for (int i = 0; i < numDevices; i++) {
      if (!AudioDevice(i).hasOutput()) continue;

      state.devices.push_back(AudioDevice(i).name());
      if (currentAudioDeviceOut == AudioDevice(i).name()) {
        state.currentDeviceOut = dev_out_index;
        state.currentOut = getLeadChannelOut() + 1;
        state.currentMaxOut = AudioDevice(i).channelsOutMax();
      }
      dev_out_index++;
    }
  };

  auto updateInDevices = [&](AudioIOState &state) {
    state.devices.clear();
    int numDevices = AudioDevice::numDevices();
    int dev_in_index = 0;
    for (int i = 0; i < numDevices; i++) {
      if (!AudioDevice(i).hasInput()) continue;

      state.devices.push_back(AudioDevice(i).name());
      if (currentAudioDeviceIn == AudioDevice(i).name()) {
        state.currentDeviceIn = dev_in_index;
        state.currentIn = getLeadChannelIn() + 1;
        state.currentMaxIn = AudioDevice(i).channelsInMax();
      }
      dev_in_index++;
    }
  };

  static std::map<AudioIO *, AudioIOState> stateMap;
  if (stateMap.find(io) == stateMap.end()) {
    stateMap[io] = AudioIOState();
    updateOutDevices(stateMap[io]);
    updateInDevices(stateMap[io]);
  }
  AudioIOState &state = stateMap[io];
  ImGui::PushID(std::to_string((unsigned long)io).c_str());

  if (io->isOpen()) {
    std::string text;
    text += "Output Device: " + state.devices.at(state.currentDeviceOut);
    text += "\nInput Device: " + state.devices.at(state.currentDeviceIn);
    text += "\nSampling Rate: " + std::to_string(int(io->fps()));
    text += "\nBuffer Size: " + std::to_string(io->framesPerBuffer());
    text += "\nOutput Channels: " + std::to_string(state.currentOut) + ", " +
            std::to_string(state.currentOut + 1);
    text += "\nInput Channels: " + std::to_string(state.currentIn) + ", " +
            std::to_string(state.currentIn + 1);
    ImGui::Text("%s", text.c_str());
    if (ImGui::Button("Stop")) {
      isPaused = true;
      io->stop();
      io->close();
      state.currentSr = getSampleRateIndex();
    }
  } else {
    if (ImGui::Button("Update Devices")) {
      updateOutDevices(state);
      updateInDevices(state);
    }

    ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
    if (ImGui::Combo("Output Device", &state.currentDeviceOut, ParameterGUI::vector_getter,
                     static_cast<void *>(&state.devices), state.devices.size())) {
      state.currentMaxOut =
        AudioDevice(state.devices.at(state.currentDeviceOut), AudioDevice::OUTPUT).channelsOutMax();
    }
    std::string chan_label_out =
      "Select Outs: (Up to " + std::to_string(state.currentMaxOut) + " )";
    ImGui::Text(chan_label_out.c_str(), "%s");
    // ImGui::SameLine();
    // ImGui::Checkbox("Mono/Stereo", &isStereo);
    // ImGui::Indent(25 * fontScale);
    // ImGui::PushItemWidth(50 * fontScale);
    ImGui::DragInt("Chan 1 out", &state.currentOut, 1.0f, 0, state.currentMaxOut - 1, "%d", 1 << 4);

    if (state.currentOut > state.currentMaxOut - 1) state.currentOut = state.currentMaxOut - 1;
    if (state.currentOut < 1) state.currentOut = 1;

    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    for (int i = 1; i < MAX_AUDIO_OUTS; i++) {
      ImGui::SameLine();
      int temp = state.currentOut + i;
      std::string channel = "Chan " + std::to_string(i + 1);
      ImGui::DragInt(channel.c_str(), &temp, 1.0f, 0, state.currentMaxOut, "%d", 1 << 4);
    }
    ImGui::PopStyleVar();

    // ImGui::Unindent(25 * fontScale);
    ImGui::PopItemWidth();

    ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x);
    if (ImGui::Combo("Input Device", &state.currentDeviceIn, ParameterGUI::vector_getter,
                     static_cast<void *>(&state.devices), state.devices.size())) {
      state.currentMaxIn =
        AudioDevice(state.devices.at(state.currentDeviceIn), AudioDevice::INPUT).channelsInMax();
    }
    std::string chan_label_in = "Select Ins: (Up to " + std::to_string(state.currentMaxIn) + " )";
    ImGui::Text(chan_label_in.c_str(), "%s");
    // ImGui::SameLine();
    // ImGui::Checkbox("Mono/Stereo", &isStereo);
    // ImGui::Indent(25 * fontScale);
    // ImGui::PushItemWidth(50 * fontScale);
    ImGui::DragInt("Chan 1 in", &state.currentIn, 1.0f, 0, state.currentMaxIn - 1, "%d", 1 << 4);

    if (state.currentIn > state.currentMaxIn - 1) state.currentIn = state.currentMaxIn - 1;
    if (state.currentIn < 1) state.currentIn = 1;

    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    for (int i = 1; i < MAX_AUDIO_INS; i++) {
      ImGui::SameLine();
      int temp = state.currentIn + i;
      std::string channel = "Chan " + std::to_string(i + 1);
      ImGui::DragInt(channel.c_str(), &temp, 1.0f, 0, state.currentMaxIn, "%d", 1 << 4);
    }
    ImGui::PopStyleVar();

    // ImGui::Unindent(25 * fontScale);
    ImGui::PopItemWidth();

    std::vector<std::string> samplingRates{"44100", "48000", "88200", "96000"};
    ImGui::Combo("Sampling Rate", &state.currentSr, ParameterGUI::vector_getter,
                 static_cast<void *>(&samplingRates), samplingRates.size());
    ImGui::PopItemWidth();
    if (ImGui::Button("Start")) {
      globalSamplingRate = std::stof(samplingRates[state.currentSr]);
      io->framesPerSecond(globalSamplingRate);
      io->framesPerBuffer(BLOCK_SIZE);
      io->deviceOut(AudioDevice(state.devices.at(state.currentDeviceOut), AudioDevice::OUTPUT));
      currentAudioDeviceOut = state.devices.at(state.currentDeviceOut);
      setOutChannels(state.currentOut - 1, state.currentMaxOut);
      io->deviceIn(AudioDevice(state.devices.at(state.currentDeviceIn), AudioDevice::INPUT));
      currentAudioDeviceIn = state.devices.at(state.currentDeviceIn);
      setInChannels(state.currentIn - 1, state.currentMaxIn);
      io->open();
      io->start();
      isPaused = false;
    }
    ImGui::SameLine();
  }
  ImGui::PopID();
}

int main() {
  AudioDevice dev = AudioDevice::defaultOutput();
  dev.print();

  CHON app;
  app.configureAudio(dev, dev.defaultSampleRate(), 1024, dev.channelsOutMax(), dev.channelsInMax());
  app.start();
  return 0;
}