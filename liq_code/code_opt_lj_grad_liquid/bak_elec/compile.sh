mpicxx -O2 -xsse3 -o grad_liquid_e_v mol-ff.cpp ff.cpp ComputePme.cpp error.cpp lattice.cpp PmeBase.cpp PmeKSpace.cpp PmeRealSpace.cpp pmetest.cpp ~/tools/fftw-2.1.5/librfftw.a ~/tools/fftw-2.1.5/libfftw.a -lm -I~/tools/fftw-2.1.5/fftw

