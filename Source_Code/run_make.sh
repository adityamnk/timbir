# !/bin/bash

#cd MBIR_4D
make clean -C ./MBIR_4D 
make -j2 -C ./MBIR_4D 
#cd ../reconstruct/basic
make clean -C ./reconstruct/basic 
make -j2 -C ./reconstruct/basic 
