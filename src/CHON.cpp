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

  execDir = al::File::directory(getExecutablePath());
  userPath = getUserHomePath();
#ifdef __APPLE__  // Uses userPath
  configFile = consts::DEFAULT_CONFIG_FILE;
  presetsPath = consts::DEFAULT_PRESETS_PATH;
  midiPresetsPath = consts::DEFAULT_MIDI_PRESETS_PATH;
  oscPresetsPath = consts::DEFAULT_OSC_PRESETS_PATH;

  samplePresetsPath = consts::DEFAULT_SAMPLE_PRESETS_PATH;

  execDir = util::getContentPath_OSX(execDir);
  al::Dir::make(userPath + consts::PERSISTENT_DATA_PATH);
  al::Dir::make(userPath + consts::DEFAULT_PRESETS_PATH);
  al::Dir::make(userPath + consts::DEFAULT_SOUND_OUTPUT_PATH);
  al::Dir::make(userPath + consts::DEFAULT_CONFIG_PATH);
  opener = "open ";
#endif

#ifdef __linux__
  std::string configPath = "/.config/CHON";
  // use xdg directories if available
  if (getenv("$XDG_CONFIG_HOME") != NULL) {
    configPath = getenv("$XDG_CONFIG_HOME") + std::string("/CHON");
  }

  configFile = configPath + "/config/config.json";
  presetsPath = configPath + "/presets";

  // create config directories if needed
  al::Dir::make(userPath + configPath + "/config");
  al::Dir::make(userPath + presetsPath);
#endif

#ifdef _WIN32
  configFile = consts::DEFAULT_CONFIG_FILE;
  presetsPath = consts::DEFAULT_PRESETS_PATH;
  midiPresetsPath = consts::DEFAULT_MIDI_PRESETS_PATH;
  oscPresetsPath = consts::DEFAULT_OSC_PRESETS_PATH;
  samplePresetsPath = consts::DEFAULT_SAMPLE_PRESETS_PATH;

  al::Dir::make(userPath + consts::PERSISTENT_DATA_PATH);
  al::Dir::make(userPath + consts::DEFAULT_PRESETS_PATH);
  al::Dir::make(userPath + consts::DEFAULT_SOUND_OUTPUT_PATH);
  al::Dir::make(userPath + consts::DEFAULT_CONFIG_PATH);
#endif

  initJsonConfig();
  json config = jsonReadConfig();
  setSoundOutputPath(config.at(consts::SOUND_OUTPUT_PATH_KEY));
  setAudioSettings(config.at(consts::SAMPLE_RATE_KEY));
  setWindowDimensions(config.at(consts::WINDOW_WIDTH_KEY), config.at(consts::WINDOW_HEIGHT_KEY));
  setFirstLaunch(config.at(consts::IS_FIRST_LAUNCH_KEY));
  setAudioDevice(config.at(consts::DEFAULT_AUDIO_DEVICE_KEY));
  setInitFullscreen(false);

  // Set output directory for presets.
  // Set output directory of recorded files.

  mPresets = std::make_unique<al::PresetHandler>(al::File::conformPathToOS(userPath + presetsPath));
  auto path = mPresets->buildMapPath("default");  // build default map if it doesn't exist
  std::ifstream exist(path);
  if (!exist.good()) {
    std::ofstream file;
    file.open(path, std::ios::out);
    file.close();
  }

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
          for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
            particleNetwork(x, y, z).resetPos(1, 0, 0);
            particleNetwork(x, y, z).resetVelocity(1, 0, 0);
          }
        }
      }
    }
  });
  yFree.registerChangeCallback([&](bool y) {
    if (!y) {
      for (int x = 1; x <= particleNetwork.sizeX(); x++) {
        for (int y = 1; y <= particleNetwork.sizeY(); y++) {
          for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
            particleNetwork(x, y, z).resetPos(0, 1, 0);
            particleNetwork(x, y, z).resetVelocity(0, 1, 0);
          }
        }
      }
    }
  });
  zFree.registerChangeCallback([&](bool z) {
    if (!z) {
      for (int x = 1; x <= particleNetwork.sizeX(); x++) {
        for (int y = 1; y <= particleNetwork.sizeY(); y++) {
          for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
            particleNetwork(x, y, z).resetPos(0, 0, 1);
            particleNetwork(x, y, z).resetVelocity(0, 0, 1);
          }
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
            for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
              particleNetwork(x, y, z).setTuningRatio(x * y);
              particleNetwork(x, y, z).setAmplitude(1.0f / x * y);
            }
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
    for (int x = 1; x <= particleNetwork.sizeX(); x++) {
      for (int y = 1; y <= particleNetwork.sizeY(); y++) {
        for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
          particleNetwork(x, y, z).setFMModFreq(newFMFreq);
        }
      }
    }
  });
  auto a_d = AudioDevice(currentAudioDeviceOut, AudioDevice::OUTPUT);
  if (!a_d.valid()) {
    audioIO().deviceOut(-1);
    currentAudioDeviceOut = AudioDevice::defaultOutput().name();
    setOutChannels(0, audioIO().channelsOutDevice());
  } else {
    audioIO().deviceOut(a_d);
    setOutChannels(config.at(consts::LEAD_CHANNEL_KEY), audioIO().channelsOutDevice());
    // TODO, make sure only 2 channels are open corresponding to out channels
    // audioIO().channelsOut({(int)config.at(consts::LEAD_CHANNEL_KEY),(int)config.at(consts::LEAD_CHANNEL_KEY)
    // + 1});
  }
  audioIO().setStreamName("CHON");
  audioIO().append(mRecorder);
  audioIO().channelsIn(0);
}

void CHON::onCreate() {  // Called when graphics context is available
  std::cout << "onCreate()" << std::endl;

  imguiInit();
  ImGuiIO &io = ImGui::GetIO();
  io.IniFilename = NULL;

  // Set if fullscreen or not.
  fullScreen(isFullScreen);

  chonReset();

  bodyFont = ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF(&RobotoMedium_compressed_data,
                                                                  RobotoMedium_compressed_size, 16);
  titleFont = ImGui::GetIO().Fonts->AddFontFromMemoryCompressedTTF(
    &RobotoMedium_compressed_data, RobotoMedium_compressed_size, 20);

  navControl().useMouse(false);
  navControl().disable();
}

void CHON::chonReset() {
  std::cout << "Reset Particle network to " << xParticles << " by " << yParticles << " by "
            << zParticles << std::endl;
  resetLock.lock();

  // Complicated logic for changing camera if switching between 1D and 2D
  if ((particleNetwork.sizeY() == 1 && yParticles > 1) && (zParticles == 1)) {
    onResize(width(), height());
  } else if ((particleNetwork.sizeZ() == 1 && zParticles > 1) && (yParticles == 1)) {
    onResize(width(), height());
  } else if (yParticles == 1 && zParticles == 1) {
    onResize(width(), height());
  }

  // reset particle network
  particleNetwork.resize(xParticles, yParticles, zParticles);

  // do I need to delete the previous chonBundles?
  xSpringGUI.reset(new ChonBundle(1));
  ySpringGUI.reset(new ChonBundle(1));
  zSpringGUI.reset(new ChonBundle(1));

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
  for (int i = 0; i <= particleNetwork.sizeZ(); i++) {
    // Create element
    auto *newSpring = new Spring{"Z Springs"};
    particleNetwork.zSprings.push_back(newSpring);
    // Register its parameter bundle with the ControlGUI
    *zSpringGUI << newSpring->bundle;
  }

  driveParticleXLeft.max(particleNetwork.sizeX());
  driveParticleYLeft.max(particleNetwork.sizeY());
  driveParticleZLeft.max(particleNetwork.sizeY());
  driveParticleXRight.max(particleNetwork.sizeX());
  driveParticleYRight.max(particleNetwork.sizeY());
  driveParticleZRight.max(particleNetwork.sizeY());

  fftDivision = log(fftBuffer.size()) /
                (particleNetwork.sizeX() * particleNetwork.sizeY() * particleNetwork.sizeZ());
  client.send("/Nparticles/",
              particleNetwork.sizeX() * particleNetwork.sizeY() * particleNetwork.sizeZ());
  resetLock.unlock();
}

void CHON::onAnimate(double dt) {  // Called once before drawing

  if (particleNetwork.sizeX() != xParticles) chonReset();

  if (particleNetwork.sizeY() != yParticles) chonReset();

  if (particleNetwork.sizeZ() != zParticles) chonReset();

  w = width();
  h = height();

  // nav().faceToward(Vec3f(0.0, 0.0, 0.0));

  particleNetwork.setFreedom(xFree, yFree, zFree);

  if (!pause) {
    if (inputOn) {
      switch (inputMode) {
        case 0:  // Peak
          driveForceLeft = inLeft.at(inLeft.getTail());
          if (driveForceLeft > inputThreshold)
            particleNetwork(driveParticleXLeft, driveParticleYLeft, driveParticleZLeft)
              .addAcceleration(driveAxisLeft, driveForceLeft * inputScale);
          if (driveStereoSplit) {
            driveForceRight = inRight.at(inRight.getTail());
            if (driveForceRight > inputThreshold)
              particleNetwork(driveParticleXRight, driveParticleYRight, driveParticleZRight)
                .addAcceleration(driveAxisRight, driveForceRight * inputScale);
          }
          break;
        case 1:  // RMS
          driveForceLeft = inLeft.getRMS(rmsSize);
          if (driveForceLeft > inputThreshold)
            particleNetwork(driveParticleXLeft, driveParticleYLeft, driveParticleZLeft)
              .addAcceleration(driveAxisLeft, driveForceLeft * inputScale);
          if (driveStereoSplit) {
            driveForceRight = inRight.getRMS(rmsSize);
            if (driveForceRight > inputThreshold)
              particleNetwork(driveParticleXRight, driveParticleYRight, driveParticleZRight)
                .addAcceleration(driveAxisRight, driveForceRight * inputScale);
          }
          break;
        case 2:  // Frequency
          particleNetwork(1, 1, 1).addAcceleration(driveAxisLeft, fftBuffer[0] * inputScale);
          for (unsigned int i = 1; i < fftBuffer.size(); i++) {
            fftIterator = floor((log(i) / fftDivision) + 1);
            particleNetwork(fftIterator % (particleNetwork.sizeX() + 1),
                            ceil(fftIterator / float(particleNetwork.sizeX() + 1)),
                            ceil(fftIterator / float(particleNetwork.sizeY() + 1)))
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

    if (dt > 0.5) {
      std::cerr << "Delta T too high for simulation" << std::endl;
    } else {
      particleNetwork.updateVelocities(dt);
    }
    for (int x = 1; x <= particleNetwork.sizeX(); x++)
      for (int y = 1; y <= particleNetwork.sizeY(); y++) {
        for (int z = 1; z <= particleNetwork.sizeZ(); z++) {  // animate stuff
          if (particleNetwork(x, y, z).particle.selected) {
            particleNetwork(x, y, z).resetVelocity(1, 1, 1);
          } else {  // increment positions according to velocity
            particleNetwork(x, y, z).velocityStep();
          }

          particleNetwork(x, y, z).updateDisplacement();

          if (DrawGraph) {  // updating the 2d displacement graph
            particleNetwork(x, y, z).graphMesh.translate(-w / (60 * graphSpeed), 0,
                                                         0);  // move previous graph data left
            if (particleNetwork(x, y, z).graphMesh.vertices().size() > w)
              particleNetwork(x, y, z).graphMesh.vertices().erase(
                particleNetwork(x, y, z).graphMesh.vertices().begin());  // erase left-edge graph
                                                                         // data
            float graphY =
              h *
              particleNetwork(x, y, z).getDisplacement()[graphAxis];  // calculate new phase value
            graphY /= particleNetwork.springLength * 8;
            graphY += (h * ((graphSpread * x / (particleNetwork.sizeX() + 1)) +
                            ((0.75 - (graphSpread / 2)))));  // add offset for each graph
            particleNetwork(x, y, z).graphMesh.vertex(w, graphY,
                                                      0);  // update with new graph data on right
          }

          if (fm) {  // copy displacement values to FM synthesis engine
            particleNetwork(x, y, z).fmSmooth.setTarget(
              particleNetwork(x, y, z).getDisplacement()[fmAxis] / particleNetwork.springLength);
          }

          if (am) {  // copy displacement values to AM synthesis engine
            float val =
              particleNetwork(x, y, z).getDisplacement()[amAxis] / particleNetwork.springLength;
            val = val > 1 ? 1 : val;
            val = val < -1 ? -1 : val;
            particleNetwork(x, y, z).amSmooth.setTarget(val);
          }

          if (stereoOn) {  // set the particle's pan position
            particleNetwork(x, y, z).panSmooth.setTarget(particleNetwork(x, y, z).x() + 0.5);
          }

          particleNetwork(x, y, z).checkZeroCrossing();

          // OSC
          if (oscOn) {
            if (oscX)
              client.send(particleNetwork(x, y, z).oscName["Xdisp"],
                          float(particleNetwork(x, y, z).getDisplacement()[0] /
                                particleNetwork.springLength));
            if (oscY)
              client.send(particleNetwork(x, y, z).oscName["Ydisp"],
                          float(particleNetwork(x, y, z).getDisplacement()[1] /
                                particleNetwork.springLength));
            if (oscZ)
              client.send(particleNetwork(x, y, z).oscName["Zdisp"],
                          float(particleNetwork(x, y, z).getDisplacement()[2] /
                                particleNetwork.springLength));
            if (oscPan)
              client.send(particleNetwork(x, y, z).oscName["Xpos"],
                          float(particleNetwork(x, y, z).x()));
          }
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
    if (particleNetwork.sizeY() > 1) {
      for (int x = 1; x <= particleNetwork.sizeX(); x++) {
        for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
          particleNetwork(x, 0, z).draw(g);
          particleNetwork(x, particleNetwork.sizeY() + 1, z).draw(g);
        }
      }
    }
    if (particleNetwork.sizeZ() > 1) {
      for (int x = 1; x <= particleNetwork.sizeX(); x++) {
        for (int y = 1; y <= particleNetwork.sizeY(); y++) {
          particleNetwork(x, y, 0).draw(g);
          particleNetwork(x, y, particleNetwork.sizeZ() + 1).draw(g);
        }
      }
    }
    if (particleNetwork.sizeX() > 1) {
      for (int y = 1; y <= particleNetwork.sizeY(); y++) {
        for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
          particleNetwork(0, y, z).draw(g);
          particleNetwork(particleNetwork.sizeX() + 1, y, z).draw(g);
        }
      }
    }
  }

  if (drawParticles) {
    for (int x = 1; x <= particleNetwork.sizeX(); x++) {
      for (int y = 1; y <= particleNetwork.sizeY(); y++) {
        for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
          g.color(
            HSV((float(x * y) / (particleNetwork.sizeX() * particleNetwork.sizeY())), 0.5, 1));
          particleNetwork(x, y, z).resetPos(!xFree, !yFree, !zFree);
          particleNetwork(x, y, z).draw(g);
        }
      }
    }
    // texBlur.copyFrameBuffer();
  }

  if (DrawGraph) {
    g.pushCamera();
    g.camera(Viewpoint::ORTHO_FOR_2D);  // Ortho [0:width] x [0:height]
    g.lighting(false);
    for (int x = 1; x <= particleNetwork.sizeX(); x++) {
      for (int y = 1; y <= particleNetwork.sizeY(); y++) {
        for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
          g.color(HSV((float(x * y * z) / (particleNetwork.sizeX() * particleNetwork.sizeY() *
                                           particleNetwork.sizeZ())),
                      0.5, 1));
          g.draw(particleNetwork(x, y, z).graphMesh);
        }
      }
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
    if (xParticles * yParticles * zParticles > 100)
      xParticles = int(100 / yParticles / zParticles);
    else if (xParticles < 1)
      xParticles = 1;
    ImGui::SameLine();
    ImGui::InputInt("y", &yParticles);
    if (xParticles * yParticles * zParticles > 100)
      yParticles = int(100 / xParticles / zParticles);
    else if (yParticles < 1)
      yParticles = 1;
    ImGui::SameLine();
    ImGui::InputInt("z", &zParticles);
    if (xParticles * yParticles * zParticles > 100)
      zParticles = int(100 / yParticles / xParticles);
    else if (zParticles < 1)
      zParticles = 1;
    ImGui::PopItemWidth();
    ImGui::Separator();
    ImGui::PushItemWidth(200);
    xSpringGUI->drawBundleGUI();
    if (particleNetwork.sizeY() > 1) ySpringGUI->drawBundleGUI();
    if (particleNetwork.sizeZ() > 1) zSpringGUI->drawBundleGUI();
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
      ImGui::Indent(20);
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
      ImGui::Indent(-20);
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
    if (ImGui::CollapsingHeader("Recorder")) {
      drawRecorderWidget(&mRecorder, audioIO().framesPerSecond(),
                         audioIO().channelsOut() <= 1 ? 1 : 2, soundOutput);
      if (ImGui::Button("Change Output Path")) {
        result = NFD_PickFolder(NULL, &outPath);

        if (result == NFD_OKAY) {
          std::string temp = outPath;
          jsonWriteToConfig(temp, consts::SOUND_OUTPUT_PATH_KEY);
          setSoundOutputPath(outPath);
        }
      }
    }
    if (ImGui::CollapsingHeader("Input")) {
      ParameterGUI::drawParameterBool(&inputOn);
      ParameterGUI::drawMenu(&inputMode);
      if (inputMode != 2) {
        ParameterGUI::drawParameterBool(&driveStereoSplit);
        if (driveStereoSplit) ImGui::Text("Left Channel");
        ParameterGUI::drawParameterInt(&driveParticleXLeft, "");
        ParameterGUI::drawParameterInt(&driveParticleYLeft, "");
        ParameterGUI::drawParameterInt(&driveParticleZLeft, "");
      }
      ParameterGUI::drawMenu(&driveAxisLeft);
      if (driveStereoSplit && inputMode != 2) {
        ImGui::Text("Right Channel");
        ParameterGUI::drawParameterInt(&driveParticleXRight, "");
        ParameterGUI::drawParameterInt(&driveParticleYRight, "");
        ParameterGUI::drawParameterInt(&driveParticleZRight, "");
        ParameterGUI::drawMenu(&driveAxisRight);
      }
      ParameterGUI::drawParameter(&inputScale);
      if (inputMode != 2) ParameterGUI::drawParameter(&inputThreshold);
      if (inputMode == 1) ParameterGUI::drawParameterInt(&rmsSize, "");
    }

    if (ImGui::CollapsingHeader("Audio IO Settings")) {
      drawAudioIO(&audioIO());
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
        "/FirstArg/SecondArg Value\n\n"
        "FirstArg options:\n"
        "pan, dispX, dispY, or dispZ \n\n"
        "SecondArg:\n"
        "Integer representing the particle number\n"
        "Value is a float\n\n"
        "Whenever the number of particles changes,\n"
        "a message is sent with the new number of particles:\n"
        "/Nparticles int");
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
    ImGui::SetNextWindowSize(ImVec2(280, 0));
    if (ImGui::BeginPopup("rightClickParticle")) {
      ImGui::PushItemWidth(ImGui::GetContentRegionAvailWidth() - 75);
      if (ImGui::InputFloat("Freq", &rightClickedFreq, 1, 10, "%.3f"))
        rightClickedParticle->setFreq(rightClickedFreq);
      if (ImGui::InputInt("Scale Step", &rightClickedStep, 1, 10)) {
        rightClickedParticle->setScaleStep(rightClickedStep);
        rightClickedParticle->setTuningRatio(particleNetwork.getScaleRatio(rightClickedStep));
        rightClickedParticle->setFreq(rightClickedParticle->getTuningRatio() *
                                      particleNetwork.getTuningRoot());
      }
      ImGui::Separator();
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
        for (int x = 1; x <= particleNetwork.sizeX(); x++) {
          for (int y = 1; y <= particleNetwork.sizeY(); y++) {
            for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
              particleNetwork(x, y, z).panSmooth.process();
            }
          }
        }
      }

      if (additiveSynthOn) {
        for (int x = 1; x <= particleNetwork.sizeX(); x++) {
          for (int y = 1; y <= particleNetwork.sizeY(); y++) {
            for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
              if (particleNetwork(x, y, z).getFreq() < 20000) {
                if (fm) {
                  particleNetwork(x, y, z).fmProcess(fmWidth);
                } else {
                  particleNetwork(x, y, z).fmProcess(0);
                }
                double sampleToAdd =
                  particleNetwork(x, y, z).processOscillator() * additiveVolume;  // scale
                sampleToAdd /= 5.0;  // scale for higher pitches
                if (am) {
                  sampleToAdd *=
                    particleNetwork(x, y, z).amSmooth.process();  // amplitude modulation
                }
                if (stereoOn) {
                  s[0] += (1 - particleNetwork(x, y, z).panSmooth.getCurrentValue()) * sampleToAdd;
                  s[1] += particleNetwork(x, y, z).panSmooth.getCurrentValue() * sampleToAdd;
                } else {
                  s[0] += sampleToAdd;
                  s[1] += sampleToAdd;
                }
              }
            }
          }
        }
      }

      if (bellSynthOn) {
        for (int x = 1; x <= particleNetwork.sizeX(); x++) {
          for (int y = 1; y <= particleNetwork.sizeY(); y++) {
            for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
              if (particleNetwork(x, y, z).isZeroTrigger(bellAxis)) {
                if (particleNetwork(x, y, z).getFreq() < 20000)
                  particleNetwork(x, y, z).bellTrigger();
              }

              double sampleToAdd = particleNetwork(x, y, z).bellProcess() * 0.2 * bellVolume;
              if (stereoOn) {
                s[0] += (1 - particleNetwork(x, y, z).panSmooth.getCurrentValue()) * sampleToAdd;
                s[1] += particleNetwork(x, y, z).panSmooth.getCurrentValue() * sampleToAdd;
              } else {
                s[0] += sampleToAdd;
                s[1] += sampleToAdd;
              }
            }
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
            for (unsigned int i = 0; i < stft.numBins(); i++) {
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
    } else
      nav().pullBack(1.2 + pow(0.002 * (1900 - width), 2));
  } else {
    nav().home();
    nav().pos(-0.5 * (500 / float(width)), 0.2, 0);
    if (yParticles > 1) {
      nav().pullBack(2.7 + pow(0.002 * (1900 - width), 2));
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
  for (int x = 1; x <= particleNetwork.sizeX(); x++) {
    for (int y = 1; y <= particleNetwork.sizeY(); y++) {
      for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
        particleNetwork(x, y, z).particle.event(PickEvent(Point, r));
      }
    }
  }
  return true;
}
bool CHON::onMouseDown(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int x = 1; x <= particleNetwork.sizeX(); x++) {
    for (int y = 1; y <= particleNetwork.sizeY(); y++) {
      for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
        if (m.button() == 0)
          particleNetwork(x, y, z).particle.event(PickEvent(Pick, r));
        else if (particleNetwork(x, y, z).particle.hover && m.button() == 2) {
          isRightClickedParticle = true;
          rightClickedParticle = &particleNetwork(x, y, z);
        }
      }
    }
  }
  return true;
}
bool CHON::onMouseDrag(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int x = 1; x <= particleNetwork.sizeX(); x++) {
    for (int y = 1; y <= particleNetwork.sizeY(); y++) {
      for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
        particleNetwork(x, y, z).particle.event(PickEvent(Drag, r));
      }
    }
  }
  return true;
}
bool CHON::onMouseUp(const Mouse &m) {
  Rayd r = getPickRay(m.x(), m.y());
  for (int x = 1; x <= particleNetwork.sizeX(); x++) {
    for (int y = 1; y <= particleNetwork.sizeY(); y++) {
      for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
        particleNetwork(x, y, z).particle.event(PickEvent(Unpick, r));
        particleNetwork(x, y, z).particle.clearSelection();
      }
    }
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
    case 'b':
      nav().print();
      return false;
    case 'g':
      drawGUI = 1 - drawGUI;
      onResize(width(), height());
      return false;
    case 'r':
      srand(std::time(0));
      for (int x = 1; x <= particleNetwork.sizeX(); x++) {
        for (int y = 1; y <= particleNetwork.sizeY(); y++) {
          for (int z = 1; z <= particleNetwork.sizeZ(); z++) {
            particleNetwork(x, y, z).addVelocity((float(rand() - (RAND_MAX / 2)) / RAND_MAX) / 10,
                                                 (float(rand() - (RAND_MAX / 2)) / RAND_MAX) / 10,
                                                 (float(rand() - (RAND_MAX / 2)) / RAND_MAX) / 10);
          }
        }
      }
      return false;
    case Keyboard::UP:
      nav().spinR(0.02);
      return false;
    case Keyboard::DOWN:
      nav().spinR(-0.02);
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
      nav().spinR(0);
      return false;
    case Keyboard::DOWN:
      nav().spinR(0);
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

void CHON::drawRecorderWidget(al::OutputRecorder *recorder, double frameRate, uint32_t numChannels,
                              std::string directory, uint32_t bufferSize) {
  struct SoundfileRecorderState {
    bool recordButton;
    bool overwriteButton;
  };
  static std::map<SoundFileBufferedRecord *, SoundfileRecorderState> stateMap;
  if (stateMap.find(recorder) == stateMap.end()) {
    stateMap[recorder] = SoundfileRecorderState{0, false};
  }
  SoundfileRecorderState &state = stateMap[recorder];
  ImGui::PushID(std::to_string((unsigned long)recorder).c_str());
  ImGui::Text("Output File Name:");
  static char buf1[64] = "test.wav";
  ImGui::PushItemWidth(ImGui::GetContentRegionAvail().x - 10.0f);
  ImGui::InputText("##Record Name", buf1, 63);
  ImGui::PopItemWidth();

  if (state.recordButton) {
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.9, 0.3, 0.3, 1.0));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.8, 0.5, 0.5, 1.0));
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0, 1.0, 1.0, 1.0));
  }
  std::string buttonText = state.recordButton ? "Stop" : "Record";
  bool recordButtonClicked = ImGui::Button(buttonText.c_str());
  if (state.recordButton) {
    ImGui::PopStyleColor();
    ImGui::PopStyleColor();
    ImGui::PopStyleColor();
  }
  if (recordButtonClicked) {
    state.recordButton = !state.recordButton;
    if (state.recordButton) {
      uint32_t ringBufferSize;
      if (bufferSize == 0) {
        ringBufferSize = 8192;
      } else {
        ringBufferSize = bufferSize * numChannels * 4;
      }
      std::string filename = buf1;
      if (!state.overwriteButton) {
        int counter = 1;
        while (File::exists(directory + filename) && counter < 9999) {
          filename = buf1;
          int lastDot = filename.find_last_of(".");
          filename = filename.substr(0, lastDot) + "_" + std::to_string(counter++) +
                     filename.substr(lastDot);
        }
      }
      if (!recorder->start(directory + filename, frameRate, numChannels, ringBufferSize,
                           gam::SoundFile::WAV, gam::SoundFile::FLOAT)) {
        std::cerr << "Error opening file for record" << std::endl;
      }
    } else {
      recorder->close();
    }
  }
  ImGui::SameLine();
  ImGui::Checkbox("Overwrite", &state.overwriteButton);
  ImGui::Text("Writing to:");
  ImGui::TextWrapped("%s", directory.c_str());

  ImGui::PopID();
}

void CHON::onExit() {
  jsonWriteToConfig(windowWidth, consts::WINDOW_WIDTH_KEY);
  jsonWriteToConfig(windowHeight, consts::WINDOW_HEIGHT_KEY);
  jsonWriteToConfig(isFullScreen, consts::FULLSCREEN_KEY);
  jsonWriteToConfig(false, consts::IS_FIRST_LAUNCH_KEY);
}

int main() {
  AudioDevice dev = AudioDevice::defaultOutput();
  dev.print();
  CHON app;
  app.start();
  return 0;
}