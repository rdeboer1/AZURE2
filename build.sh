#!/bin/bash


# remove build directory
rm -r build

#check if the build directory exists or not
#if not, then create the build directory
[ ! -d /build ] && mkdir -p ./build

#go into the build directory
cd ./build

#run cmake and if sucessfull run make install 
cmake .. -DUSE_QWT=ON && make && make install
