#!/bin/bash

git submodule update --init --recursive

(
  mkdir -p build
  cd build
  mkdir -p release
  cd release
  cmake -DCMAKE_BUILD_TYPE=Release -Wno-deprecated -DBUILD_EXAMPLES=0 ../..
)

# Configure debug build
(
  mkdir -p build
  cd build
  mkdir -p debug
  cd debug
  cmake -DCMAKE_BUILD_TYPE=Debug -Wno-deprecated -DBUILD_EXAMPLES=0 ../..
)
