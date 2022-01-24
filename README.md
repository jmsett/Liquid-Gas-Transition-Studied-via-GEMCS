# Liquid-Gas-Transition-Studied-via-GEMCS
The main idea of the project is to investigate the liquid-gas transition of a fore-shifted Lenard-Jones particle system via Gibbs ensemble Monte Carlo simulations.
In particular, we wish to determine the system's liquid-gas binodal in the temperature-density plane and the phase diagram in the pressure-temperature plane.
We will compare the obtained binodal with results from the literature, and we will discuss the behavior around the critical point in light of the van der Waals
equation of state.

The full report can be found in the repository as a PDF file. The plots can be found in the 'GEMCS plots.ipynb' jupyter-notebook (new path names are needed in order to load the data).

## Program
### main.cpp:
Main GEMCS file.
The following simulation parameters can be changed here: totV, V_BoxA, V_boxB, Ntot, N_BoxA, N_BoxB, temperature, dmax, Vmax, eqsweeps (number of sweeps for thermalization), sweeps (GEMCS sweeps).

### global_random.h:
Global random functions.

### progresBar.h:
Progress Bar animation for the terminal output.

### particles.h:
Particle and Box properties.

### potential.h:
Functions for computing the total potential energy, as well as rescaled potential energy.

### montecarlo.h:
Monte Carlo Moves. MC move probabilities can be changed Pn, Pd, Pv.

### static_correlations.h:
Static correlation function g(r) (radial distribution function).

### verletlist3D:
Verlet list definition for faster GEMC simulations, as opposed to updating pair interactions every round. (Uses celllist.h)

### celllist.h:
Celllist for faster GEMC simulations (used by verletlist3D.h)

## Data
The folder '100K' contains the results of the GEMCS for different temperatures, 'gr' contains the results of the radial distribution function for different temperatures, and 'thermal_eq' contains the thermal equilibrium data.

### binodal.txt:
Binodal data GEMCS numerical results.

### delta_rho.dat:
rho_l - rho_g data points for fit via the scaling law.

### paper_exp.dat, paper_theo.dat:
Data points of the binodals in Verlet's paper, as mentioned in the report (experimental/numerical binodals).

### pressure.dat:
Pressure numerical data for the gas and liquid phase.

### p_T.dat:
Pressure data of the GEMCS for different temperatures, taking the average of liquid-gas pressure numerical results from pressure.dat.

### rho_avg.dat:
data points for the liquid-gas average density for fit with the coexistence width line.

### shifted.dat:
Data points of the binodal in Vrabec's paper.

### sigma_binodal.dat:
Data on the Gaussian fit's standard deviation for the points on the T-rho phase diagram.

### time_equilibrium.txt:
Data on the total time of the simulation for every temperature, and data on the equilibrium point for every temperature.
