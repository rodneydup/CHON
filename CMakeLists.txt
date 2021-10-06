cmake_minimum_required(VERSION 3.0)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.10" CACHE STRING "Minimum OS X deployment version")

set(APP_NAME CHON)
project(${APP_NAME} VERSION 1.2
                    DESCRIPTION "Coupled Harmonic Oscillator Network"
                    LANGUAGES CXX)

# add allolib as a subdirectory to the project
add_subdirectory(external/allolib)
add_subdirectory(external/al_ext)

if (NOT TARGET al_soundfile)
  message(FATAL_ERROR "al_soundfile not built. Likely libsndfile not found")
endif()

# Hacky (and fragile way) to move libsndfile-1.dll to bin/Resources
if(WIN32)
  get_filename_component(sndfile_dir ${SNDFILE_LIBRARY} DIRECTORY)
  set(${sndfile_dir}/../bin/libsndfile-1.dll)
  message("WINDOWS: Attempting to copy sndfile dll to bin/")
  file(COPY ${sndfile_dir}/../bin/libsndfile-1.dll
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)
endif()

# OSX library needed for native file dialog.
if(APPLE)
  find_library(AppKit Appkit)
endif()

# Linux library needed for native file dialog.
if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(GTK3 REQUIRED gtk+-3.0)
  include_directories(${GTK3_INCLUDE_DIRS})
  link_directories(${GTK3_LIBRARY_DIRS})
endif()

set(includes
src/structs.hpp
src/CHON.hpp
src/Roboto-Medium.hpp
src/consts.hpp
)

set(srcs
src/CHON.cpp
)


set(LIB_NAME chonLib)
add_library(${LIB_NAME}
  STATIC 
    ${includes}
    ${srcs}
)


target_link_libraries(${LIB_NAME} PUBLIC 
  al 
)

target_include_directories(${LIB_NAME} PRIVATE 
  ${CMAKE_SOURCE_DIR}/external
)

set_target_properties(${LIB_NAME} PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
)

# bunch of stuff for packaging MacOS Apps
if(APPLE)
  find_library(AppKit Appkit)
  target_link_libraries(${LIB_NAME} PUBLIC ${Appkit})

  set(RESOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin/CHON.app/Contents/Resources")
  set(CMAKE_MACOSX_RPATH 1)

  set(MACOSX_BUNDLE_ICON_FILE "CHON.icns")
  set(ICON ${CMAKE_CURRENT_SOURCE_DIR}/deployment/icons/CHON.icns)
  set_source_files_properties(${ICON} PROPERTIES
       MACOSX_PACKAGE_LOCATION "Resources")

  set(CMAKE_INSTALL_RPATH "@executable_path/../Resources/")

  set(RESOURCE_LIST
    "libsndfile"
  )

  foreach(RESOURCE_DIR_NAME ${RESOURCE_LIST})
    SET(RESOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deployment/externalResources/${RESOURCE_DIR_NAME}")
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      file(GLOB CHILDREN RELATIVE "${RESOURCE_DIR}" ${RESOURCE_DIR}/*)
      foreach(CHILD ${CHILDREN})
        IF(IS_DIRECTORY "${RESOURCE_DIR}/${CHILD}")
          # be wary of .DS_Store files. glob them out by using [A-Za-z]*.*
          FILE(GLOB RESOURCE_DIR_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${RESOURCE_DIR}/${CHILD}/[A-Za-z]*.*")
          SET_SOURCE_FILES_PROPERTIES(
            ${RESOURCE_DIR_FILES}
            PROPERTIES MACOSX_PACKAGE_LOCATION Resources/${RESOURCE_DIR_NAME}/${CHILD}
          )
          SET(RESOURCE_FILES ${RESOURCE_FILES} "deployment/externalResources/${RESOURCE_DIR_FILES}")
        ELSEIF(NOT ${CHILD} MATCHES "\\.DS_Store|.*\\.am")
          SET_SOURCE_FILES_PROPERTIES(
            "${RESOURCE_DIR}/${CHILD}"
            PROPERTIES MACOSX_PACKAGE_LOCATION Resources/${RESOURCE_DIR_NAME}
          )
          SET(RESOURCE_FILES ${RESOURCE_FILES} "deployment/externalResources/${RESOURCE_DIR_NAME}/${CHILD}")
        ENDIF()
    endforeach()
  else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") #WORK ON
    INSTALL( DIRECTORY ${RESOURCE_DIR} DESTINATION "${DATA_DIR}" PATTERN "*.am" EXCLUDE)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  endforeach()

  add_executable(${APP_NAME} MACOSX_BUNDLE
    src/CHON.cpp
    ${RESOURCE_FILES}
    ${ICON}
  )
  
endif(APPLE)

if(WIN32)
  set(ICON ${CMAKE_CURRENT_SOURCE_DIR}/deployment/icons/CHON-windows.rc)
  # path to main source file
  add_executable(${APP_NAME} 
    src/CHON.cpp
    ${ICON}
  )
endif(WIN32)

if(UNIX AND NOT APPLE)
  # path to main source file
  add_executable(${APP_NAME} src/CHON.cpp)
endif()

target_include_directories(${APP_NAME} PRIVATE
${CMAKE_SOURCE_DIR}/external
${CMAKE_SOURCE_DIR}/external/nativefiledialog/src/include
)

find_library(NFD_LIBRARY nfd HINTS 
  ${CMAKE_SOURCE_DIR}/external/nativefiledialog/build/lib/Release/x64
)

# if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/al_ext)
#   message("Buiding extensions in al_ext")
#   add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/al_ext)
#   get_target_property(AL_EXT_LIBRARIES al_ext AL_EXT_LIBRARIES)
#   target_link_libraries(${APP_NAME} PRIVATE ${AL_EXT_LIBRARIES})
# endif()

# link allolib to project
get_target_property(AL_EXT_LIBRARIES al_ext AL_EXT_LIBRARIES)

target_link_libraries(${APP_NAME} PUBLIC 
  al 
  "${AL_EXT_LIBRARIES}"
  ${NFD_LIBRARY}
) 
if(APPLE)
  target_link_libraries(${APP_NAME} PUBLIC ${Appkit})
endif()

if(LINUX)
  target_link_libraries(${APP_NAME} PUBLIC ${GTK3_LIBRARIES})
endif()

# Link windows libraries.
if(WIN32)
  target_link_libraries(${APP_NAME} PUBLIC comctl32)
endif()


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

# Link windows libraries.
if(WIN32)
  target_link_libraries(${LIB_NAME} PUBLIC comctl32)
endif()

target_link_libraries(${APP_NAME} PUBLIC ${LIB_NAME} ${LIBRARIES})

# binaries are put into the ./bin directory by default
set_target_properties(${APP_NAME} PROPERTIES
  BUILD_WITH_INSTALL_RPATH 1
  INSTALL_NAME_DIR "@rpath"
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
  add_custom_command (TARGET ${APP_NAME}
    POST_BUILD 
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libsndfile/lib/libsndfile.1.dylib" 
      "@rpath/libsndfile/libsndfile.1.dylib"
      "$<TARGET_FILE:${APP_NAME}>" 

    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libvorbis/lib/libvorbis.0.dylib" 
      "@rpath/libsndfile/libvorbis.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libsndfile.1.dylib"
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/flac/lib/libFLAC.8.dylib" 
      "@rpath/libsndfile/libFLAC.8.dylib"
      "${RESOURCE_PATH}/libsndfile/libsndfile.1.dylib"
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libogg/lib/libogg.0.dylib" 
      "@rpath/libsndfile/libogg.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libsndfile.1.dylib"
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libvorbis/lib/libvorbisenc.2.dylib" 
      "@rpath/libsndfile/libvorbisenc.2.dylib"
      "${RESOURCE_PATH}/libsndfile/libsndfile.1.dylib" 

    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
    "-change" "/usr/local/opt/libogg/lib/libogg.0.dylib" 
    "@rpath/libsndfile/libogg.0.dylib"
    "${RESOURCE_PATH}/libsndfile/libFLAC.8.dylib"

    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libogg/lib/libogg.0.dylib" 
      "@rpath/libsndfile/libogg.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libvorbis.0.dylib"

    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/opt/libogg/lib/libogg.0.dylib" 
      "@rpath/libsndfile/libogg.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libvorbisenc.2.dylib" 
    COMMAND "${CMAKE_INSTALL_NAME_TOOL}"
      "-change" "/usr/local/Cellar/libvorbis/1.3.6/lib/libvorbis.0.dylib" 
      "@rpath/libsndfile/libvorbis.0.dylib"
      "${RESOURCE_PATH}/libsndfile/libvorbisenc.2.dylib" 
    VERBATIM
  )

endif(APPLE)
