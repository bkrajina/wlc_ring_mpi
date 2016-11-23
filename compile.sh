#/bin/bash
cd src
mpifort -o wlcsim  mersenne_twister.f90 SIMcode/inputparams.f90 SIMcode/params.f03 SIMcode/main.f03 SIMcode/*.f95 MCcode/*.f95 /usr/lib/liblapack.so -O5
rm ../data/*
rm ../data/*/*
cd ../
mv src/wlcsim ./

