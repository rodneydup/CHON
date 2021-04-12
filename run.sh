#!/bin/bash

THREADS=$1
if [[ "$THREADS" == "" ]]; then
  echo "
   NOTE: you can append '-j#', replacing # with your CPU's number of cores +1, to compile faster.
   
   Example: './scripts/run.sh -j5' (for a 4-core processor)
   "
fi

(
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cmake --build ./build/release $THREADS
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    cmake --build ./build/release $THREADS
  elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    echo "BUILD WINDOWS TODO"
    cmake --build ./build/release $THREADS
  fi
)

result=$?
if [ ${result} == 0 ]; then
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    ./bin/CHON
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    open bin/CHON.app
  elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    echo "RUN WINDOWS TODO"
    ./EmissionControl2.exe -DRTAUDIO_API_JACK=1 -DRTMIDI_API_JACK=0
  fi
fi
