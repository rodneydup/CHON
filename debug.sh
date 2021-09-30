#!/bin/bash

THREADS=$1
if [[ "$THREADS" == "" ]]; then
  echo "
   NOTE: you can append '-j#', replacing # with your CPU's number of cores +1, to compile faster.
   
   Example: './scripts/debug.sh -j5' (for a 4-core processor)
   "
fi

(
  if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    cmake --build ./build/debug --config Debug $THREADS
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    cmake --build ./build/debug --config Debug $THREADS
  elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    cmake --build ./build/debug --config Debug $THREADS
  fi
)
result=$?
if [ ${result} == 0 ]; then
  if [ $(uname -s) == "Linux" ]; then
    gdb -ex run ./CHON
  elif [ $(uname -s) == "Darwin" ]; then
    chmod 444 ./CHON.app/Contents/Resources/libsndfile/*
    lldb ./CHON.app/Contents/MacOS/CHON
  elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    echo "RUN WINDOWS DEBUGGER TODO"
  fi

fi
