# NVT &amp; NPT Monte Carlo simulation of LJ gas

### NVT ensemble:

Compile as: `make lj-nvt`

To run: `./lj-nvt.x [Char(20) UUID of the run] [Int #Atoms] [Real Temp[Kelvin]] [Real SimBoxSize[Angstrom]] [Int #EQsteps] [Int #PRODsteps] [Int samplingRate] [Int seed]`

Example:
`./lj-nvt.x 'testRun-' 100 300.0 10.0 100000 100000 100 3127387 > sim.out`

Starts simulation of 100 particles interacting through LJ potential at temperature 300 Kelvin in square box with volume 10.0 x 10.0 x 10.0 Angstrom^3. The simulation
first performs 10^5 equilibration sweeps, each consisting of 100 local moves follwed by 10^5 production sweeps during which the observables are sampled every 100 production sweeps. The sampled values of observables (Energy, Pressure, Virial) are written to `sim.out`. Samples of configuration obtained during production sweeps are stored in file `testRun-conf.xyz` in .xyz file format with distances in LJ units.

### NPT ensemble:

Compile as: `make lj-npt`

To run: `./lj-npt-cr.x [Char(20) UUID of the run] [Int #Atoms] [Real Temp[Kelvin]] [Real Pressure[LJ Units]] [Int #EQsteps] [Int #PRODsteps] [Int samplingRate] [Int seed]`

Example:
`./lj-npt-cr.x 'testRun-' 100 300.0 1.8 100000 100000 100 3127387 > sim.out`

Starts simulation of 100 particles interacting through LJ potential at temperature 300 Kelvin at reduced pressure 1.8 (In LJ Units). The simulation
first performs 10^5 equilibration sweeps, each consisting of 100 local moves 
and single volume scaling move follwed by 10^5 production sweeps during which the observables are sampled every 100 production sweeps. The sampled values of observables (Energy, Density, Virial) are written to `sim.out`. Samples of configuration obtained during production sweeps are stored in file `testRun-conf.xyz` in .xyz file format with distances in LJ units.

### Data analysis:

##### Sample means

`bining-n[v|p]t.py` computes averages and errors of sampled observables from NVT or NPT simulation by decorrelating the data through binning analysis. 

To run: `python3 binning-n[v|p]t.py [String filename] [Int MaxBinSize]`

Example: `python3 binning-nvt.py 'sim.out' 50 > averages.dat`

##### Radial distribution
Compile as: `make rdist`

To run: `./rDist.x [Char(20) filename] [Int #Samples] [Real maxR[LJ units]] [Int #Bins]`

Example: `./rDist.x testRun-conf.xyz 1000 2.5 60`

Computes the radial distribution function from 1000 samples of configuration stored in file `testRun-conf.xyz` with pair separation cutoff set to 2.5 in LJ units with resolution of 60 bins. The result is stored in `rho_r.dat`.