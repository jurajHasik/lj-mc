lj-nvt:
	rm -f lj-nvt.x
	gfortran -fno-range-check -O5 -o lj-nvt.x MersenneTwister.f90 lj-NVT.f90
	rm -f *.mod
lj-nvt-RLX77:
	rm -f lj-nvt-RLX77.x
	gfortran -O3 -c ranlux_cernlib.F
	gfortran -O3 -o lj-nvt-RLX77.x ranlux_cernlib.o lj-NVT-RLX77.f90
	rm -f *.mod
lj-nvt-RLX90:
	rm -f lj-nvt-RLX90.x
	gfortran -O3 -o lj-nvt-RLX90.x ranlux.f90 lj-NVT-RLX90.f90
	rm -f *.mod
lj-npt:
	gfortran -fno-range-check -O5 -o lj-npt-cr.x MersenneTwister.f90 lj-NPT-cutRscl.f90
	gfortran -fno-range-check -O5 -o lj-npt-ncr.x MersenneTwister.f90 lj-NPT-noCutRscl.f90
	rm -f *.mod
rDist:
	rm -f rDist.x
	gfortran -fno-range-check -o rDist.x radial-dist.f90
	rm -f *.mod
clean:
	rm -f *.mod *.x
