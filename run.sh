#!/bin/bash

THREADS=$1
if [[ "$THREADS" == "" ]]; then
  echo "
   NOTE: you can append '-j#', replacing # with your CPU's number of cores +1, to compile faster.
   
   Example: './scripts/run.sh -j5' (for a 4-core processor)
   "
fi

(
  # utilizing cmake's parallel build options
  # recommended: -j <number of processor cores + 1>
  # This is supported in cmake >= 3.12 use -- -j5 for older versions
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
  ./bin/CHON
fi
