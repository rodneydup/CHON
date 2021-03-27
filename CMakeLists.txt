cmake_minimum_required(VERSION 3.0)

project(CHON)

# name of application. replace 'app' with desired app name
set(APP_NAME CHON)

# bunch of stuff for packaging MacOS Apps
if(APPLE)
  set(RESOURCE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin/CHON.app/Contents/Resources")
  set(CMAKE_MACOSX_RPATH 1)

  set(MACOSX_BUNDLE_ICON_FILE "CHON.icns")
  set(ICON ${CMAKE_CURRENT_SOURCE_DIR}/deployment/icons/CHON.icns)
  set_source_files_properties(${ICON} PROPERTIES
       MACOSX_PACKAGE_LOCATION "Resources")

  set(CMAKE_INSTALL_RPATH "@executable_path/../Resources/")
  set(RESOURCE_LIST
  )

  foreach(RESOURCE_DIR_NAME ${RESOURCE_LIST})
    SET(RESOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../externalResources/${RESOURCE_DIR_NAME}")
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
          SET(RESOURCE_FILES ${RESOURCE_FILES} "../externalResources/${RESOURCE_DIR_FILES}")
        ELSEIF(NOT ${CHILD} MATCHES "\\.DS_Store|.*\\.am")
          SET_SOURCE_FILES_PROPERTIES(
            "${RESOURCE_DIR}/${CHILD}"
            PROPERTIES MACOSX_PACKAGE_LOCATION Resources/${RESOURCE_DIR_NAME}
          )
          SET(RESOURCE_FILES ${RESOURCE_FILES} "../externalResources/${RESOURCE_DIR_NAME}/${CHILD}")
        ENDIF()
      endforeach()
    else(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") #WORK ON
      INSTALL( DIRECTORY ${RESOURCE_DIR} DESTINATION "${DATA_DIR}" PATTERN "*.am" EXCLUDE)
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  endforeach()
endif(APPLE)

# path to main source file
add_executable(${APP_NAME} src/CHON.cpp)

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
