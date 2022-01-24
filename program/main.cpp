#include <iostream>
#include <random>
#include "static_correlations.h"
#include "particles.h"
#include "potential.h"
#include "montecarlo.h"
#include "global_random.h"
#include "verletlist3D.h"
#include <chrono>

int main(int argc, char **argv)
{
  // Define total (fixed) Volume and Box volume fraction.
  double totV = 1500.0; // 0.5 density ->1500
  double V_BoxA = totV/2.0;
  double V_BoxB = totV - V_BoxA;
  // Define total (fixed) particle-number and Box particle-number.
  int Ntot = 750;
  int N_BoxA = Ntot/2.0;
  int N_BoxB = Ntot - N_BoxA;

  // Define temperature.
  double temperature = 0.8;
  // Define Canonical MC displacement parameter "dmax".
  double dmax = 0.05;
  // Define Volume Exchange Parameter "Vmax".
  double Vmax = 100.0; // try 0.1 computing p

  // Define Particles in Box A and Box B.
  Particles BoxA(N_BoxA, V_BoxA);
  Particles BoxB(N_BoxB, V_BoxB);
  //Particles BoxA("eq_Box_V750_N375.dat");
  //Particles BoxB("eq_Box_V750_N375.dat");

  // Start at random positions.
  BoxA.setRandomPositions();
  BoxB.setRandomPositions();

  // Define Potential
  Potential::WCA_U potential(1.0,2.5);

  // Define VerletList3D for both Boxes
  VerletList3D vlistA(BoxA.n_particles_, BoxA.x_, BoxA.y_, BoxA.z_, BoxA.lx_, potential.cutoff_);
  VerletList3D vlistB(BoxB.n_particles_, BoxB.x_, BoxB.y_, BoxB.z_, BoxB.lx_, potential.cutoff_);

  std::cout << "T = " << temperature << std::endl;
  std::cout << "[ MC Simulation: Started ]" << std::endl;

  // Reach equilibrium with Canonical MC
  int eqsweeps = 4000; //10000
  //BoxA.writeState("iState.dat");
  double canRate;
  Montecarlo::Canonical(BoxA, potential, eqsweeps*(BoxA.n_particles_*1.0/(Ntot*1.0)), dmax, temperature, vlistA,canRate);
  //BoxA.writeState("eq_Box_V750_N375.dat");
  //BoxA.writeState("fState.dat");

  // Calculate the correlation function gr for Box A
  CorrelationsFunctions::calculate_gr(200, 1000, BoxA, potential, dmax, temperature, vlistA);

  Montecarlo::Canonical(BoxB, potential, eqsweeps*(BoxB.n_particles_*1.0/(Ntot*1.0)), dmax, temperature, vlistB,canRate);

  auto start = std::chrono::high_resolution_clock::now();

  int sweeps = 200000; //1500
  double d_accRate = 0.0;
  double v_accRate = 0.0;
  double p_accRate = 0.0;
  Montecarlo::Gibbs(BoxA, BoxB, potential, sweeps, dmax, Vmax, temperature, vlistA, vlistB, d_accRate, v_accRate, p_accRate);

  std::cout << "[ MC Simulation: Finished ]" << std::endl;
  std::cout << "Rates:  d: " << d_accRate << ", v: " << v_accRate << ", p: " << p_accRate << std::endl;

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count() << " s\n";
  return 0;
}
