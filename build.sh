mkdir build && cd build 
cmake ../AZURE2/ -DBUILD_GUI=ON -DUSE_QWT=ON && make -j4 && cd -

mkdir -p $PREFIX/bin
cp build/src/AZURE2 $PREFIX/bin/