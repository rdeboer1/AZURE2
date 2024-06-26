#!/bin/bash


# remove build directory
rm -r build

#check if the build directory exists or not
#if not, then create the build directory
[ ! -d /build ] && mkdir -p ./build

#go into the build directory
cd ./build

#run cmake and if sucessfull run make install 
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_QWT=ON -DBUILD_LIBRARY=OFF -DCMAKE_INSTALL_PREFIX=/opt/Azure2 -DMINUIT2_INCLUDE_DIR=/usr/include/Minuit2 -DMINUIT2_LIBRARY_DIR=/usr/lib/root -DQWT_VERSION_STRING=6.2.0 -DQt5_DIR=/usr/lib/cmake/Qt5 -DQt5Widgets_DIR=/usr/lib/cmake/Qt5Widgets -DQt5PrintSupport_DIR=/usr/lib/cmake/Qt5PrintSupport -DQt5Gui_DIR=/usr/lib/cmake/Qt5Gui -DQt5Core_DIR=/usr/lib/cmake/Qt5Core -DROOT_DIR=/usr/lib/cmake/ROOT -Dnlohmann_json_DIR=/usr/share/cmake/nlohmann_json -DBUILD_LIBRARY=OFF "$PWD/.."
make
