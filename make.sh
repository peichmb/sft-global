module load intel
rm *.o
rm *.mod
clear
ifort -O3 -c -traceback -g -pg -check bounds iobelda.F90\
      params.F90 random.F90 grid.F90 lsflows.F90 vars.F90 output.F90 \
      sources.F90 derivatives.F90 rk4.F90 timestep.F90 crank.F90 inflows.F90 \
      initconds.F90 update.F90 cpoint.F90 init.F90 \
      xerbla.f zgtsv.f dgtsv.f randgen.f -lfftw3
ifort -O3 -traceback -g -pg -check bounds \
      params.o random.o cpoint.o grid.o lsflows.o vars.o sources.o derivatives.o \
      initconds.o rk4.o timestep.o crank.o inflows.o update.o iobelda.o \
      output.o init.o xerbla.o zgtsv.o dgtsv.o randgen.o main.F90 -lfftw3 -o sft.x
echo 'Done'
