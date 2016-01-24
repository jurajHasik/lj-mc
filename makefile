lj-nvt:
	rm -f lj-nvt.x
	gfortran -fno-range-check -o lj-nvt.x MersenneTwister.f90 lj-NVT.f90
	rm -f *.mod
lj-npt:
	rm -f lj-npt.x
	gfortran -fno-range-check -o lj-npt.x MersenneTwister.f90 lj-NPT.f90
	rm -f *.mod
rDist:
	rm -f rDist.x
	gfortran -fno-range-check -o rDist.x radial-dist.f90
	rm -f *.mod
clean:
	rm -f *.mod *.x
