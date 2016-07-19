FC=gfortran

MT_FLAGS=-fno-range-check

MT=prng/MersenneTwister.f90
# Suppres warning during compilation of MersenneTwister.f90
# See https://gcc.gnu.org/onlinedocs/gfortran/Fortran-Dialect-Options.html
MT_FLAGS=-fno-range-check

RLX90=prng/ranlux.f90
RLX77=prng/ranlux_cernlib.F

lj-nvt: clean
	rm -f lj-nvt.x
	$(FC) $(MT_FLAGS) -O3 -o lj-nvt.x $(MT) lj-NVT.f90
	rm -f *.mod
lj-nvt-RLX77: clean
	rm -f lj-nvt-RLX77.x
	$(FC) -O3 -c $(RLX77)
	$(FC) -O3 -o lj-nvt-RLX77.x ranlux_cernlib.o lj-NVT-RLX77.f90
	rm -f *.mod
lj-nvt-RLX90: clean
	rm -f lj-nvt-RLX90.x
	$(FC) -O3 -o lj-nvt-RLX90.x $(RLX90) lj-NVT-RLX90.f90
	rm -f *.mod
lj-npt: clean
	$(FC) $(MT_FLAGS) -O3 -o lj-npt-cr.x $(MT) lj-NPT-cutRscl.f90
	# Unstable wrt. volume scaling move
	#$(FC) $(MT_FLAGS) -O3 -o lj-npt-ncr.x $(MT) lj-NPT-noCutRscl.f90
	rm -f *.mod
rDist: clean
	rm -f rDist.x
	$(FC) -o rDist.x data-analysis/radial-dist.f90
	rm -f *.mod
clean:
	rm -f *.mod *.x *.o
