#!/bin/bash

#check if the build directory exists or not
#if not, then create the build directory
[ ! -d /build ] && mkdir -p ./build

#go into the build directory
cd ./build

#run cmake and if sucessfull run make install 
cmake .. && make && make install
