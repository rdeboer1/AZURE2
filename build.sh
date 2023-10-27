mkdir build && cd build 
cmake .. -DBUILD_GUI=ON -DUSE_QWT=ON && make -j4 && cd -

mkdir -p $PREFIX/bin
cp build/src/AZURE2 $PREFIX/bin/
