/bin/rm *.o
gfortran -c -O3 -mcmodel=medium func.f pip.f pythag.f svbksb.f svdcmp.f svdfit.f
gfortran *.o -o pip.x
