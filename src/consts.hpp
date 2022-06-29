#ifndef CONSTS_HPP
#define CONSTS_HPP

#include <string>

namespace consts {

const std::string VERSION_NUMBER = "1.2";

const std::string SOUND_OUTPUT_PATH_KEY = "USER_SOUND_OUTPUT_PATH";
const std::string SAMPLE_RATE_KEY = "SAMPLE_RATE";
const std::string BUFFER_SIZE_KEY = "BUFFER_SIZE";
const std::string WINDOW_WIDTH_KEY = "WINDOW_WIDTH";
const std::string WINDOW_HEIGHT_KEY = "WINDOW_HEIGHT";
const std::string IS_FIRST_LAUNCH_KEY = "FIRST_LAUNCH";
const std::string DEFAULT_AUDIO_DEVICE_OUT_KEY = "DEFAULT_AUDIO_DEVICE_OUT";
const std::string DEFAULT_AUDIO_DEVICE_IN_KEY = "DEFAULT_AUDIO_DEVICE_IN";
const std::string LEAD_CHANNEL_OUT_KEY = "LEAD_CHANNEL_OUT";
const std::string LEAD_CHANNEL_IN_KEY = "LEAD_CHANNEL_IN";

const int SAMPLE_RATE = 48000;
const int BUFFER_SIZE = 1024;
const int MAX_AUDIO_OUTS = 2;
const int MAX_AUDIO_INS = 2;
const float WINDOW_WIDTH = 1920;
const float WINDOW_HEIGHT = 1080;
const bool IS_FIRST_LAUNCH = true;
const std::string DEFAULT_AUDIO_DEVICE_OUT = "";
const std::string DEFAULT_AUDIO_DEVICE_IN = "";
const size_t DEFAULT_LEAD_CHANNEL_OUT = 0;
const size_t DEFAULT_LEAD_CHANNEL_IN = 0;

/**
 *  DEFAULT USER PATHS
 */

#ifdef __APPLE__
const std::string PERSISTENT_DATA_PATH = "/Music/CHON/";
const std::string DEFAULT_SOUND_OUTPUT_PATH = PERSISTENT_DATA_PATH + "soundOutput/";
const std::string DEFAULT_PRESETS_PATH = PERSISTENT_DATA_PATH + "presets/";
const std::string DEFAULT_CONFIG_PATH = PERSISTENT_DATA_PATH + "configs/";
const std::string DEFAULT_CONFIG_FILE = DEFAULT_CONFIG_PATH + "config.json";

#endif

#ifdef __linux__
const std::string DEFAULT_SOUND_OUTPUT_PATH = "";
#endif

#ifdef _WIN32
const std::string PERSISTENT_DATA_PATH = "/Music/CHON/";
const std::string DEFAULT_SOUND_OUTPUT_PATH = PERSISTENT_DATA_PATH + "soundOutput/";
const std::string DEFAULT_PRESETS_PATH = PERSISTENT_DATA_PATH + "presets/";
const std::string DEFAULT_CONFIG_PATH = PERSISTENT_DATA_PATH + "configs/";
const std::string DEFAULT_CONFIG_FILE = DEFAULT_CONFIG_PATH + "config.json";
#endif

}  // namespace consts
#endif