#!/bin/bash

git submodule update --init --recursive

rm -f external/al_ext/assets3d/CMakeLists.txt
rm -f external/al_ext/openvr/CMakeLists.txt
rm -f external/al_ext/spatialaudio/CMakeLists.txt
rm -f external/al_ext/statedistribution/CMakeLists.txt

# Build nativefiledialog if it doesnt exist../external/libsamplerate/build
(
  if [ ! -d "./external/nativefiledialog/build/lib" ]; then
    cd external/nativefiledialog/build/
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
      cd gmake_linux
      make config=release_x64
    elif [[ "$OSTYPE" == "darwin"* ]]; then # note: can't get make file to work, relies on xcode bleh
      cd gmake_macosx
      export MACOSX_DEPLOYMENT_TARGET=10.10
      make config=release_x64
      #cd xcode4
      #xcodebuild -scheme nfd build -project nfd.xcodeproj/ -configuration Release CFLAGS=-mmacosx-version-min=10.10 CXXFLAGS=-mmacosx-version-min=10.10
    elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
      #DEPENDENT ON VISUAL STUDIO
      cd vs2010/
      msbuild.exe ./nfd.vcxproj //p:Configuration=Release //p:Platform=x64
    fi
  fi
)

JACKSUPPORT=1
# Check flags to see if user wants to remove JACK support for pipewire compatibility
PIPEWIRE=$1

if [[ "$PIPEWIRE" == "" ]]; then
  echo "

  CONFIGURING EC2 WITH JACK SUPPORT. THERE MAY BE ISSUES ON SYSTEMS RUNNING PIPEWIRE!
  To build without jack for pipewire compatibility, run run configure with the argument "pipewire" (i.e. ./configure.sh pipewire)

  "
elif [[ "$PIPEWIRE" == "pipewire" ]]; then
  echo "

  CONFIGURING EC2 WITHOUT JACK FOR PIPEWIRE COMPATIBILITY.

  "
  JACKSUPPORT=0
fi

# CMake configure Release

(
  mkdir -p build
  cd build
  mkdir -p release
  cd release
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cmake -DCMAKE_BUILD_TYPE=Release -Wno-deprecated -DBUILD_EXAMPLES=0 -DRTAUDIO_API_JACK=$JACKSUPPORT -DRTMIDI_API_JACK=1 ../..
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
