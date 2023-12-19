mkdir build && cd build 
cmake .. && make -j6 && cd -

mkdir -p $PREFIX/bin
cp build/src/AZURE2 $PREFIX/bin/
