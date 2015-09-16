# !/bin/bash

cd MBIR_4D
make clean
make -j2
cd ../reconstruct/basic
make clean 
make -j2