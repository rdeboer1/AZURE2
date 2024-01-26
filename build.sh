mkdir build && cd build && cmake -DCMAKE_PREFIX_PATH=$CONDA_PREFIX .. && make -j6 && cd -

mkdir -p $PREFIX/bin
cp build/src/AZURE2 $PREFIX/bin/