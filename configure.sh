#!/bin/bash

git submodule update --init --recursive

rm -f external/al_ext/assets3d/CMakeLists.txt
rm -f external/al_ext/openvr/CMakeLists.txt
rm -f external/al_ext/spatialaudio/CMakeLists.txt
rm -f external/al_ext/statedistribution/CMakeLists.txt

(
  mkdir -p build
  cd build
  mkdir -p release
  cd release
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cmake -DCMAKE_BUILD_TYPE=Release -Wno-deprecated -DBUILD_EXAMPLES=0 -DRTAUDIO_API_JACK=1 -DRTMIDI_API_JACK=1 ../..
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    cmake -DCMAKE_BUILD_TYPE=Release -Wno-deprecated -DBUILD_EXAMPLES=0 -DRTAUDIO_API_JACK=0 -DRTMIDI_API_JACK=0 ../..
  elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    cmake -DCMAKE_BUILD_TYPE=Release -Wno-deprecated -DBUILD_EXAMPLES=0 -DRTAUDIO_API_JACK=0 -DRTMIDI_API_JACK=0 ../..
  fi
)

# Configure debug build
(
  cd build
  mkdir -p debug
  cd debug
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cmake -DCMAKE_BUILD_TYPE=Debug -Wno-deprecated -DBUILD_EXAMPLES=0 -DRTAUDIO_API_JACK=1 -DRTMIDI_API_JACK=1 ../..
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    cmake -DCMAKE_BUILD_TYPE=Debug -Wno-deprecated -DBUILD_EXAMPLES=0 -DRTAUDIO_API_JACK=0 -DRTMIDI_API_JACK=0 ../..
  elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    cmake -DCMAKE_BUILD_TYPE=Debug -Wno-deprecated -DBUILD_EXAMPLES=0 -DRTAUDIO_API_JACK=0 -DRTMIDI_API_JACK=0 ../..
  fi
)
