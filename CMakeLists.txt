cmake_minimum_required(VERSION 3.0)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.10" CACHE STRING "Minimum OS X deployment version")

set(APP_NAME CHON)
project(${APP_NAME} VERSION 1.0
                    DESCRIPTION "Coupled oscillator signal generator"
                    LANGUAGES CXX)



# bunch of stuff for packaging MacOS Apps
if(APPLE)
  find_library(AppKit Appkit)

  set(RESOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin/CHON.app/Contents/Resources")
  set(CMAKE_MACOSX_RPATH 1)

  set(MACOSX_BUNDLE_ICON_FILE "CHON.icns")
  set(ICON ${CMAKE_CURRENT_SOURCE_DIR}/deployment/icons/CHON.icns)
  set_source_files_properties(${ICON} PROPERTIES
       MACOSX_PACKAGE_LOCATION "Resources")

  set(CMAKE_INSTALL_RPATH "@executable_path/../Resources/")
  add_executable(${APP_NAME} MACOSX_BUNDLE
  src/CHON.cpp
  ${RESOURCE_FILES}
  ${ICON}
  )
endif(APPLE)

if(UNIX AND NOT APPLE)
  # path to main source file
  add_executable(${APP_NAME} src/CHON.cpp)
endif()

# add allolib as a subdirectory to the project
add_subdirectory(allolib)

if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/al_ext)
  message("Buiding extensions in al_ext")
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/al_ext)
  get_target_property(AL_EXT_LIBRARIES al_ext AL_EXT_LIBRARIES)
  target_link_libraries(${APP_NAME} PRIVATE ${AL_EXT_LIBRARIES})
endif()

# link allolib to project
target_link_libraries(${APP_NAME} PRIVATE al)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # On macOS if jack is present, RtMidi will link to jack but not provide the
    # /usr/local/lib location as link directory. This statement only serves the
    # purpose of working around that bug. Perhaps can be removed on future
    # updates of RtMidi
    link_directories("/usr/local/lib")
    link_directories("/Library/Frameworks/Jackmp.framework/Versions/A/")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# example line for find_package usage
# find_package(Qt5Core REQUIRED CONFIG PATHS "C:/Qt/5.12.0/msvc2017_64/lib" NO_DEFAULT_PATH)

# replace ${PATH_TO_INCLUDE_DIR} before including other libraries
# target_include_directories(${APP_NAME} PRIVATE ${PATH_TO_INCLUDE_DIR})

# replace ${PATH_TO_LIB_FILE} before linking other libraries
# target_link_libraries(${APP_NAME} PRIVATE ${PATH_TO_LIB_FILE})

# binaries are put into the ./bin directory by default
set_target_properties(${APP_NAME} PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_LIST_DIR}/bin
  RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_LIST_DIR}/bin
)

if(APPLE)
  set_target_properties(${APP_NAME} PROPERTIES
  MACOSX_BUNDLE_BUNDLE_NAME "${APP_NAME}"
  MACOSX_BUNDLE_GUI_IDENTIFIER "${APP_NAME}"
  MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/deployment/MacOS/MacOSBundleInfo.plist.in
  )

endif(APPLE)
