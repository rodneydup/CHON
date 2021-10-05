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
      for (int x = 1; x <= particleNetwork.sizeX(); x++) {
        for (int y = 1; y <= particleNetwork.sizeY(); y++) {
          particleNetwork(x, y).resetPos(1, 0, 0);
          particleNetwork(x, y).resetVelocity(1, 0, 0);
        }
      }
    }
  });
  yFree.registerChangeCallback([&](bool y) {
    if (!y) {
      for (int x = 1; x <= particleNetwork.sizeX(); x++) {
        for (int y = 1; y <= particleNetwork.sizeY(); y++) {
          particleNetwork(x, y).resetPos(0, 1, 0);
          particleNetwork(x, y).resetVelocity(0, 1, 0);
        }
      }
    }
  });
  zFree.registerChangeCallback([&](bool z) {
    if (!z) {
      for (int x = 1; x <= particleNetwork.sizeX(); x++) {
        for (int y = 1; y <= particleNetwork.sizeY(); y++) {
          particleNetwork(x, y).resetPos(0, 0, 1);
          particleNetwork(x, y).resetVelocity(0, 0, 1);
        }
      }
    }
  });
  reverbTail.registerChangeCallback([&](float tail) { reverb.decay(tail); });

  graphAxis.setElements({"x", "y", "z"});
  fmAxis.setElements({"x", "y", "z"});
  amAxis.setElements({"x", "y", "z"});
  bellAxis.setElements({"x", "y", "z"});
  tuningScale.setElements({"Pentatonic", "Major", "Chromatic", "Harmonics", "Bohlen-Pierce"});
  driveAxisLeft.setElements({"x", "y", "z"});
  driveAxisRight.setElements({"x", "y", "z"});
  inputMode.setElements({"Peak", "RMS", "Frequency"});

  tuningScale.registerChangeCallback([&](int tuning) {
    switch (tuning) {
      case 0:
        particleNetwork.retune(pentScale);
        break;
      case 1:
        particleNetwork.retune(majScale);
        break;
      case 2:
        particleNetwork.retune(chromScale);
        break;
      case 3:
        for (int x = 1; x <= particleNetwork.sizeX(); x++) {
          for (int y = 1; y <= particleNetwork.sizeY(); y++) {
            particleNetwork(x, y).setTuningRatio(x * y);
            particleNetwork(x, y).setAmplitude(1.0f / x * y);
          }
        }
        particleNetwork.retune();
        break;
      case 4:
        particleNetwork.retune(bpScale);
        break;
      default:
        break;
    }
  });

  tuningRoot.registerChangeCallback([&](float root) {
    particleNetwork.setTuningRoot(root);
    particleNetwork.retune();
  });

  tuningScale.set(0);
  tuningRoot.set(440);

  fmFreq.registerChangeCallback([&](float newFMFreq) {
    for (int x = 1; x <= particleNetwork.sizeX(); x++)
      for (int y = 1; y <= particleNetwork.sizeY(); y++)
        particleNetwork(x, y).setFMModFreq(newFMFreq);
  });

  audioIO().setStreamName("CHON");
}

void CHON::onCreate() {  // Called when graphics context is available
  std::cout << "onCreate()" << std::endl;

  imguiInit();
  ImGuiIO &io = ImGui::GetIO();
  io.IniFilename = NULL;

  chonReset();

  bodyFont = ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF(&RobotoMedium_compressed_data,
                                                                  RobotoMedium_compressed_size, 16);
  titleFont = ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF(
    &RobotoMedium_compressed_data, RobotoMedium_compressed_size, 20);

  navControl().useMouse(false);
  navControl().disable();
}

void CHON::chonReset() {
  std::cout << "Reset Particle network to " << xParticles << " by " << yParticles << std::endl;
  resetLock.lock();

  // changing camera depending on if it's a 2D or 1D particle system
  if (particleNetwork.sizeY() == 1 && yParticles > 1) {
    onResize(width(), height());
  } else if (yParticles == 1) {
    onResize(width(), height());
  }

  // reset particle network
  particleNetwork.resize(xParticles, yParticles);

  // do I need to delete the previous chonBundles?
  xSpringGUI.reset(new ChonBundle(1));
  ySpringGUI.reset(new ChonBundle(1));

  for (int i = 0; i <= particleNetwork.sizeX(); i++) {
    // Create element
    auto *newSpring = new Spring{"X Springs"};
    particleNetwork.xSprings.push_back(newSpring);
    // Register its parameter bundle with the ControlGUI
    *xSpringGUI << newSpring->bundle;
  }
  for (int i = 0; i <= particleNetwork.sizeY(); i++) {
    // Create element
    auto *newSpring = new Spring{"Y Springs"};
    particleNetwork.ySprings.push_back(newSpring);
    // Register its parameter bundle with the ControlGUI
    *ySpringGUI << newSpring->bundle;
  }

  driveParticleXLeft.max(particleNetwork.sizeX());
  driveParticleYLeft.max(particleNetwork.sizeY());
  driveParticleXRight.max(particleNetwork.sizeX());
  driveParticleYRight.max(particleNetwork.sizeY());

  fftDivision = log(fftBuffer.size()) / (particleNetwork.sizeX() * particleNetwork.sizeY());
  client.send("/Nparticles/", particleNetwork.sizeX() * particleNetwork.sizeY());
  resetLock.unlock();
}

void CHON::onAnimate(double dt) {  // Called once before drawing

  if (particleNetwork.sizeX() != xParticles) chonReset();

  if (particleNetwork.sizeY() != yParticles) chonReset();

  w = width();
  h = height();

  particleNetwork.setFreedom(xFree, yFree, zFree);

  if (!pause) {
    if (inputOn) {
      switch (inputMode) {
        case 0:  // Peak
          driveForceLeft = inLeft.at(inLeft.getTail());
          if (driveForceLeft > inputThreshold)
            particleNetwork(driveParticleXLeft, driveParticleYLeft)
              .addAcceleration(driveAxisLeft, driveForceLeft * inputScale);
          if (driveStereoSplit) {
            driveForceRight = inRight.at(inRight.getTail());
            if (driveForceRight > inputThreshold)
              particleNetwork(driveParticleXRight, driveParticleYRight)
                .addAcceleration(driveAxisRight, driveForceRight * inputScale);
          }
          break;
        case 1:  // RMS
          driveForceLeft = inLeft.getRMS(rmsSize);
          if (driveForceLeft > inputThreshold)
            particleNetwork(driveParticleXLeft, driveParticleYLeft)
              .addAcceleration(driveAxisLeft, driveForceLeft * inputScale);
          if (driveStereoSplit) {
            driveForceRight = inRight.getRMS(rmsSize);
            if (driveForceRight > inputThreshold)
              particleNetwork(driveParticleXRight, driveParticleYRight)
                .addAcceleration(driveAxisRight, driveForceRight * inputScale);
          }
          break;
        case 2:  // Frequency
          particleNetwork(1, 1).addAcceleration(driveAxisLeft, fftBuffer[0] * inputScale);
          for (int i = 1; i < fftBuffer.size(); i++) {
            fftIterator = floor((log(i) / fftDivision) + 1);
            particleNetwork(fftIterator % (particleNetwork.sizeX() + 1),
                            ceil(fftIterator / float(particleNetwork.sizeX() + 1)))
              .addAcceleration(driveAxisLeft, fftBuffer[i] * inputScale);
          }
          fftIterator = 0;
          break;
        default:
          break;
      }
    }

    particleNetwork.setMass(mass);
    particleNetwork.setDamping(damping);

    particleNetwork.updateVelocities(dt);

    for (int x = 1; x <= particleNetwork.sizeX(); x++)
      for (int y = 1; y <= particleNetwork.sizeY(); y++) {  // animate stuff
        if (particleNetwork(x, y).particle.selected) {
          particleNetwork(x, y).resetVelocity(1, 1, 1);
        } else {  // increment positions according to velocity
          particleNetwork(x, y).velocityStep();
        }

        particleNetwork(x, y).updateDisplacement();

        if (DrawGraph) {  // updating the 2d displacement graph
          particleNetwork(x, y).graphMesh.translate(-w / (60 * graphSpeed), 0,
                                                    0);  // move previous graph data left
          if (particleNetwork(x, y).graphMesh.vertices().size() > w)
            particleNetwork(x, y).graphMesh.vertices().erase(
              particleNetwork(x, y).graphMesh.vertices().begin());  // erase left-edge graph data
          float graphY =
            h * particleNetwork(x, y).getDisplacement()[graphAxis];  // calculate new phase value
          graphY /= particleNetwork.springLength * 8;
          graphY += (h * ((graphSpread * x / (particleNetwork.sizeX() + 1)) +
                          ((0.75 - (graphSpread / 2)))));  // add offset for each graph
          particleNetwork(x, y).graphMesh.vertex(w, graphY,
                                                 0);  // update with new graph data on right
        }

        if (fm) {  // copy displacement values to FM synthesis engine
          particleNetwork(x, y).fmSmooth.setTarget(particleNetwork(x, y).getDisplacement()[fmAxis] /
                                                   particleNetwork.springLength);
        }

        if (am) {  // copy displacement values to AM synthesis engine
          float val =
            particleNetwork(x, y).getDisplacement()[amAxis] / particleNetwork.springLength;
          val = val > 1 ? 1 : val;
          val = val < -1 ? -1 : val;
          particleNetwork(x, y).amSmooth.setTarget(val);
        }

        if (stereoOn) {  // set the particle's pan position
          particleNetwork(x, y).panSmooth.setTarget(particleNetwork(x, y).x() + 0.5);
        }

        particleNetwork(x, y).checkZeroCrossing();

        // OSC
        if (oscOn) {
          if (oscX)
            client.send(
              particleNetwork(x, y).oscName["Xdisp"],
              float(particleNetwork(x, y).getDisplacement()[0] / particleNetwork.springLength));
          if (oscY)
            client.send(
              particleNetwork(x, y).oscName["Ydisp"],
              float(particleNetwork(x, y).getDisplacement()[1] / particleNetwork.springLength));
          if (oscZ)
            client.send(
              particleNetwork(x, y).oscName["Zdisp"],
              float(particleNetwork(x, y).getDisplacement()[2] / particleNetwork.springLength));
          if (oscPan)
            client.send(particleNetwork(x, y).oscName["Xpos"], float(particleNetwork(x, y).x()));
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
    if (particleNetwork.sizeY() > 1)
      for (int x = 1; x <= particleNetwork.sizeX(); x++) {
        particleNetwork(x, 0).draw(g);
        particleNetwork(x, particleNetwork.sizeY() + 1).draw(g);
      }
    for (int y = 1; y <= particleNetwork.sizeY(); y++) {
      particleNetwork(0, y).draw(g);
      particleNetwork(particleNetwork.sizeX() + 1, y).draw(g);
    }
  }

  if (drawParticles) {
    for (int x = 1; x <= particleNetwork.sizeX(); x++)
      for (int y = 1; y <= particleNetwork.sizeY(); y++) {
        g.color(HSV((float(x * y) / (particleNetwork.sizeX() * particleNetwork.sizeY())), 0.5, 1));
        particleNetwork(x, y).resetPos(!xFree, !yFree, !zFree);
        particleNetwork(x, y).draw(g);
      }
    // texBlur.copyFrameBuffer();
  }

  if (DrawGraph) {
    g.pushCamera();
    g.camera(Viewpoint::ORTHO_FOR_2D);  // Ortho [0:width] x [0:height]
    g.lighting(false);
    for (int x = 1; x <= particleNetwork.sizeX(); x++)
      for (int y = 1; y <= particleNetwork.sizeY(); y++) {
        g.color(HSV((float(x * y) / (particleNetwork.sizeX() * particleNetwork.sizeY())), 0.5, 1));
        g.draw(particleNetwork(x, y).graphMesh);
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
    if (particleNetwork.sizeY() > 1) ySpringGUI->drawBundleGUI();
    ImGui::PopItemWidth();
    ParameterGUI::drawParameter(&mass);
    ParameterGUI::drawParameter(&damping);
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
    ParameterGUI::drawMenu(&tuningScale);
    ParameterGUI::drawParameter(&tuningRoot);
    if (ImGui::CollapsingHeader("Bell Synth")) {
      ParameterGUI::drawParameterBool(&bellSynthOn);
      ImGui::SetNextItemWidth(35);
      ImGui::SameLine();
      ParameterGUI::drawMenu(&bellAxis);
      ParameterGUI::drawParameter(&bellVolume);
    }
    if (ImGui::CollapsingHeader("Additive Synth")) {
      ParameterGUI::drawParameterBool(&additiveSynthOn);
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

    if (isRightClickedParticle) {
      ImGui::OpenPopup("rightClickParticle");
      rightClickedFreq = rightClickedParticle->getFreq();
      rightClickedAmplitude = rightClickedParticle->getAmplitude();
      rightClickedStep = rightClickedParticle->getScaleStep();
      isRightClickedParticle = false;
    }
    ImGui::SetNextWindowSize(ImVec2(250, 0));
    if (ImGui::BeginPopup("rightClickParticle")) {
      ImGui::PushItemWidth(ImGui::GetContentRegionAvailWidth() - 35);
      if (ImGui::InputFloat("Freq", &rightClickedFreq, 1, 10, "%.3f"))
        rightClickedParticle->setFreq(rightClickedFreq);
      if (ImGui::InputInt("Scale Step", &rightClickedStep, 1, 10)) {
        rightClickedParticle->setScaleStep(rightClickedStep);
        rightClickedParticle->setTuningRatio(particleNetwork.getScaleRatio(rightClickedStep));
        rightClickedParticle->setFreq(rightClickedParticle->getTuningRatio() *
                                      particleNetwork.getTuningRoot());
      }
      if (ImGui::SliderFloat("Vol", &rightClickedAmplitude, 0, 1, "%.3f"))
        rightClickedParticle->setAmplitude(rightClickedAmplitude);
      ImGui::PopItemWidth();
      ImGui::EndPopup();
    }

    imguiEndFrame();

    imguiDraw();
  }
}

void CHON::onSound(AudioIOData &io) {  // Audio callback
  while (io()) {
    double s[2] = {0, 0};

    if (resetLock.try_lock()) {
      if (stereoOn) {
        for (int x = 1; x <= particleNetwork.sizeX(); x++)
          for (int y = 1; y <= particleNetwork.sizeY(); y++) {
            particleNetwork(x, y).panSmooth.process();
          }
      }

      if (additiveSynthOn) {
        for (int x = 1; x <= particleNetwork.sizeX(); x++)
          for (int y = 1; y <= particleNetwork.sizeY(); y++) {
            if (particleNetwork(x, y).getFreq() < 20000) {
              if (fm) {
                particleNetwork(x, y).fmProcess(fmWidth);
              } else {
                particleNetwork(x, y).fmProcess(0);
              }
              double sampleToAdd =
                particleNetwork(x, y).processOscillator() * additiveVolume;  // scale
              sampleToAdd /= 5.0;  // scale for higher pitches
              if (am) {
                sampleToAdd *= particleNetwork(x, y).amSmooth.process();  // amplitude modulation
              }
              if (stereoOn) {
                s[0] += (1 - particleNetwork(x, y).panSmooth.getCurrentValue()) * sampleToAdd;
                s[1] += particleNetwork(x, y).panSmooth.getCurrentValue() * sampleToAdd;
              } else {
                s[0] += sampleToAdd;
                s[1] += sampleToAdd;
              }
            }
          }
      }

      if (bellSynthOn) {
        for (int x = 1; x <= particleNetwork.sizeX(); x++)
          for (int y = 1; y <= particleNetwork.sizeY(); y++) {
            if (particleNetwork(x, y).isZeroTrigger(bellAxis)) {
              if (particleNetwork(x, y).getFreq() < 20000) particleNetwork(x, y).bellTrigger();
            }

            double sampleToAdd = particleNetwork(x, y).bellProcess() * 0.2 * bellVolume;
            if (stereoOn) {
              s[0] += (1 - particleNetwork(x, y).panSmooth.getCurrentValue()) * sampleToAdd;
              s[1] += particleNetwork(x, y).panSmooth.getCurrentValue() * sampleToAdd;
            } else {
              s[0] += sampleToAdd;
              s[1] += sampleToAdd;
            }
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
  for (int x = 1; x <= particleNetwork.sizeX(); x++)
    for (int y = 1; y <= particleNetwork.sizeY(); y++)
      particleNetwork(x, y).particle.event(PickEvent(Point, r));
  return true;
}
bool CHON::onMouseDown(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int x = 1; x <= particleNetwork.sizeX(); x++)
    for (int y = 1; y <= particleNetwork.sizeY(); y++) {
      if (m.button() == 0)
        particleNetwork(x, y).particle.event(PickEvent(Pick, r));
      else if (particleNetwork(x, y).particle.hover && m.button() == 2) {
        isRightClickedParticle = true;
        rightClickedParticle = &particleNetwork(x, y);
      }
    }
  return true;
}
bool CHON::onMouseDrag(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int x = 1; x <= particleNetwork.sizeX(); x++)
    for (int y = 1; y <= particleNetwork.sizeY(); y++)
      particleNetwork(x, y).particle.event(PickEvent(Drag, r));
  return true;
}
bool CHON::onMouseUp(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int x = 1; x <= particleNetwork.sizeX(); x++)
    for (int y = 1; y <= particleNetwork.sizeY(); y++) {
      particleNetwork(x, y).particle.event(PickEvent(Unpick, r));
      particleNetwork(x, y).particle.clearSelection();
    }
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
      for (int x = 1; x <= particleNetwork.sizeX(); x++)
        for (int y = 1; y <= particleNetwork.sizeY(); y++) {
          particleNetwork(x, y).addVelocity((float(rand() - (RAND_MAX / 2)) / RAND_MAX) / 10,
                                            (float(rand() - (RAND_MAX / 2)) / RAND_MAX) / 10,
                                            (float(rand() - (RAND_MAX / 2)) / RAND_MAX) / 10);
        }
      return false;
    case Keyboard::UP:
      nav().moveF(0.02);
      return false;
    case Keyboard::DOWN:
      nav().moveF(-0.02);
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
      nav().moveF(0);
      return false;
    case Keyboard::DOWN:
      nav().moveF(0);
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