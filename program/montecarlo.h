#pragma once

#include "particles.h"
#include "potential.h"
#include "global_random.h"
#include "verletlist3D.h"
#include <random>
#include <vector>
#include <iostream>
#include "progresBar.h"

namespace Montecarlo
{
 // Namespace function declaration
 void SingleDisplace(Particles&, Potential::WCA_U&, double, double, VerletList3D&, double&);
 void BoxSingleDisplace(Particles&, Particles&, Potential::WCA_U&, double, double, VerletList3D&, VerletList3D&, double&);
 void VolumeExchange(Particles&, Particles&, Potential::WCA_U&, double, double, VerletList3D&, VerletList3D&, double&);

 void ParticleExchange(Particles&, Particles&, Potential::WCA_U&, double, VerletList3D&, VerletList3D&, double&);
 void exchange(Particles&, Particles&, Potential::WCA_U&, double, VerletList3D&, VerletList3D&, double&);



void Gibbs(Particles &BoxA, Particles &BoxB, Potential::WCA_U &potential, int n_sweeps, double dmax, double Vmax, double temperature, VerletList3D &vlistA, VerletList3D &vlistB, double &d_accRate, double &v_accRate, double &p_accRate)
{
  progresBar bar(n_sweeps * (BoxA.n_particles_ + BoxB.n_particles_), "Gibbs MC ", "");
  int c1 = 0;
  int c2 = 0;
  int c3 = 0;

  std::ofstream outfile;
  std::ofstream epotfile;
  outfile.open("V_N.dat");
  epotfile.open("epot.dat");

  double Switch;
  double n_d = 1000.;//1
  double n_ve = 1.;//1
  double n_pe = 100.;//200, 350

  for(int j=0; j < n_sweeps; j++)
  {
    epotfile << j << " " << potential.compute_Utot(BoxA,vlistA) << " " << potential.compute_Utot(BoxB,vlistB) << " " << potential.compute_P(temperature,BoxA,vlistA) << " " << potential.compute_P(temperature,BoxB,vlistB) << std::endl;
    outfile << j << " " << BoxA.BoxV_ << " " << BoxA.n_particles_ << " " << BoxB.BoxV_ << " " << BoxB.n_particles_ << std::endl;

    for(int i=0; i < (BoxA.n_particles_ + BoxB.n_particles_); i++)
    {
      //outfile << i+j*(BoxA.n_particles_ + BoxB.n_particles_) << " " << BoxA.BoxV_ << " " << BoxA.n_particles_
      //        << " " << BoxB.BoxV_ << " " << BoxB.n_particles_ << std::endl;

      Switch = globRNG::uniform(globRNG::gen) * (n_d + n_ve + n_pe);
      if(Switch < n_d)
      {
        bar.subname("Displace");
        bar.update(c1,c2,c3);
        // Choose a random Box and displace a random particle from that Box.
        BoxSingleDisplace(BoxA, BoxB, potential, dmax, temperature, vlistA, vlistB, d_accRate);
        c1++;
      }
      else if(Switch < n_d + n_ve)
      {
        bar.subname("  Volume Exchange");
        bar.update(c1,c2,c3);
        // Randomly exchange volume dv in [-Vmax,Vmax], BoxA -> BoxB or BoxB -> BoxA.
        VolumeExchange(BoxA, BoxB, potential, Vmax, temperature, vlistA, vlistB, v_accRate);
        c2++;
      }
      else
      {
        bar.subname("Particle Exchange");
        bar.update(c1,c2,c3);
          // choose random box, eliminate random particle in I, create random particle in II.
          ParticleExchange(BoxA, BoxB, potential, temperature, vlistA, vlistB, p_accRate);
        c3++;
      }
    } // end of particle loop
  } // end of sweep loop
  d_accRate /= 1.0*c1;
  v_accRate /= 1.0*c2;
  p_accRate /= 1.0*c3;
  bar.end();
  outfile.close();
}

void ParticleExchange(Particles &BoxA, Particles &BoxB, Potential::WCA_U &potential, double temperature, VerletList3D &vlistA, VerletList3D &vlistB, double &p_accRate)
{
  // Select a random Box.
  double boxRand = globRNG::uniform(globRNG::gen);

  if(boxRand < 0.5) // Remove particle from Box A and add to Box B.
  {
    exchange(BoxA, BoxB, potential, temperature, vlistA, vlistB, p_accRate);
  }

  else // Remove particle from Box B and add to Box A.
  {
    exchange(BoxB, BoxA, potential, temperature, vlistB, vlistA, p_accRate);
  }
}

void exchange(Particles &BoxA, Particles &BoxB, Potential::WCA_U &potential, double temperature, VerletList3D &vlistA, VerletList3D &vlistB, double &p_accRate)
// Exchange Particles from BoxA to BoxB.
{
  double beta = 1.0 / temperature;
  int nA = 1.0*BoxA.n_particles_;
  int nB = 1.0*BoxB.n_particles_;

  // Select a random particle from BoxA and compute delta_epotA.
  int i = globRNG::uniform(globRNG::gen) * BoxA.n_particles_;
  double delta_epotA = - potential.compute_Ui(BoxA, i, vlistA);

  // Get random coordinates and compute virtual delta_epotB.
  double x_Rand = globRNG::uniform(globRNG::gen) * BoxB.lx_;
  double y_Rand = globRNG::uniform(globRNG::gen) * BoxB.lx_;
  double z_Rand = globRNG::uniform(globRNG::gen) * BoxB.lx_;

  double delta_epotB = potential.computeVirtual(BoxB, x_Rand, y_Rand, z_Rand);

  double Pacc = (nA/(nB + 1.0)) * (BoxB.BoxV_/BoxA.BoxV_) * exp(-beta * (delta_epotA + delta_epotB));
  double rcomp = globRNG::uniform(globRNG::gen);

  if ( (Pacc > 1) or (Pacc < 1 and rcomp < Pacc) ) // Accepted.
  {
    BoxA.removePart(i);
    BoxB.addParticle(x_Rand, y_Rand, z_Rand);
    vlistA.reset();
    vlistB.reset();
    vlistA.redefine(BoxA.n_particles_, BoxA.x_, BoxA.y_, BoxA.z_, BoxA.lx_, potential.cutoff_);
    vlistB.redefine(BoxB.n_particles_, BoxB.x_, BoxB.y_, BoxB.z_, BoxB.lx_, potential.cutoff_);
    //vlistA.n_update(BoxA.n_particles_, BoxA.x_, BoxA.y_, BoxA.z_);
    //vlistB.n_update(BoxB.n_particles_, BoxB.x_, BoxB.y_, BoxB.z_);
    p_accRate++;
  }

}

void VolumeExchange(Particles &BoxA, Particles &BoxB, Potential::WCA_U &potential, double Vmax, double temperature, VerletList3D &vlistA, VerletList3D &vlistB, double &v_accRate)
{
  double beta = 1.0 / temperature;

  // Draw random double in [-Vmax,Vmax).
  double dV = (2 * globRNG::uniform(globRNG::gen) - 1) * Vmax;

  // Save old variables
  double A_l_old = BoxA.lx_;
  double B_l_old = BoxB.lx_;
  double Box_A_V_old = BoxA.BoxV_;
  double Box_B_V_old = BoxB.BoxV_;

  // Define new Volumes, lengths.
  BoxA.BoxV_ += dV;
  BoxB.BoxV_ -= dV;
  double A_l_new = pow(BoxA.BoxV_,1.0/3.0);
  double B_l_new = pow(BoxB.BoxV_,1.0/3.0);

  // Compute delta_epot.
  double delta_epotA = -potential.compute_Utot(BoxA, vlistA);
  delta_epotA += potential.scaled_Utot(A_l_old, A_l_new);

  double delta_epotB = -potential.compute_Utot(BoxB, vlistB);
  delta_epotB += potential.scaled_Utot(B_l_old, B_l_new);

  // Define acceptance probability.
  double Pacc = pow(BoxA.BoxV_/Box_A_V_old, BoxA.n_particles_) * pow(BoxB.BoxV_/Box_B_V_old, BoxB.n_particles_);
  Pacc *= exp(-beta * (delta_epotA + delta_epotB));

  // Draw a random number in (0,1]
  double rcomp = globRNG::uniform(globRNG::gen);
  if (Pacc < 1 and rcomp > Pacc) // Rejected.
  {
    BoxA.BoxV_ = Box_A_V_old;
    BoxB.BoxV_ = Box_B_V_old;
  }
  else // Accepted.
  {
    v_accRate++;

    BoxA.l_  = A_l_new;
    BoxA.lx_ = A_l_new;
    BoxA.ly_ = A_l_new;
    BoxA.lz_ = A_l_new;

    BoxB.l_  = B_l_new;
    BoxB.lx_ = B_l_new;
    BoxB.ly_ = B_l_new;
    BoxB.lz_ = B_l_new;

    // Coordinate Scaling
    for(int i = 0; i < BoxA.n_particles_; ++i)
    {
      BoxA.x_[i] *= (A_l_new/A_l_old);
      BoxA.y_[i] *= (A_l_new/A_l_old);
      BoxA.z_[i] *= (A_l_new/A_l_old);
    }

    for(int i = 0; i < BoxB.n_particles_; ++i)
    {
      BoxB.x_[i] *= (B_l_new/B_l_old);
      BoxB.y_[i] *= (B_l_new/B_l_old);
      BoxB.z_[i] *= (B_l_new/B_l_old);
    }
    //vlistA.v_update(BoxA.lx_, BoxA.x_, BoxA.y_, BoxA.z_);
    //vlistB.v_update(BoxB.lx_, BoxB.x_, BoxB.y_, BoxB.z_);
    vlistA.reset();
    vlistB.reset();
    vlistA.redefine(BoxA.n_particles_, BoxA.x_, BoxA.y_, BoxA.z_, A_l_new, potential.cutoff_);
    vlistB.redefine(BoxB.n_particles_, BoxB.x_, BoxB.y_, BoxB.z_, B_l_new, potential.cutoff_);
  }

}

void Canonical(Particles &particles, Potential::WCA_U &potential, int n_sweeps, double dmax, double temperature, VerletList3D &vlist, double &accRate)
{
  progresBar bar(n_sweeps * particles.n_particles_, "Canonical");

  double epot;
  std::ofstream outepot;
  outepot.open("epot_can.dat");

  //assume kb = 1;
  double beta = 1.0 / temperature;

  for (int i = 0; i < n_sweeps; i++)
  {
    epot = potential.compute_Utot(particles,vlist);
    outepot << i << " " << epot << std::endl;

    for (int j = 0; j < particles.n_particles_; j++)
    {
      bar.update();
      //choose a random particle
      int p = static_cast<int>(globRNG::uniform(globRNG::gen)*particles.n_particles_);
      //compute the old energy
      double delta_epot = -potential.compute_Ui(particles, p, vlist);

      //shift the particle position
      double sx = (2 * globRNG::uniform(globRNG::gen) - 1) * dmax;
      double sy = (2 * globRNG::uniform(globRNG::gen) - 1) * dmax;
      double sz = (2 * globRNG::uniform(globRNG::gen) - 1) * dmax;
      //std::cout << sx << " " << sy<< " " << sz << std::endl;

      vlist.move(particles, p, sx, sy, sz);
      //particles.applyBoundCond(p);

      // Compute epot_new - epot_old.
      delta_epot += potential.compute_Ui(particles, p, vlist);

      //draw a random number in (0,1]
      double rcomp = globRNG::uniform(globRNG::gen);

      accRate++;

      if (rcomp > exp(-beta * delta_epot) and delta_epot > 0) //rejection
      {
        vlist.move(particles, p, -sx, -sy, -sz);
        accRate--;
      }
    }//end loop over particles
  }//end loop over n_sweeps
  accRate /= (1.0 * particles.n_particles_ * n_sweeps);
  bar.end();
  outepot.close();
}

// Canonical with no file write
void Canonical2(Particles &particles, Potential::WCA_U &potential, int n_sweeps, double dmax, double temperature, VerletList3D &vlist)
{
  //assume kb = 1;
  double beta = 1.0 / temperature;

  for (int i = 0; i < n_sweeps; i++)
  {
    for (int j = 0; j < particles.n_particles_; j++)
    {
      //choose a random particle
      int p = static_cast<int>(globRNG::uniform(globRNG::gen)*particles.n_particles_);
      //compute the old energy
      double delta_epot = -potential.compute_Ui(particles, p, vlist);

      //shift the particle position
      double sx = (2 * globRNG::uniform(globRNG::gen) - 1) * dmax;
      double sy = (2 * globRNG::uniform(globRNG::gen) - 1) * dmax;
      double sz = (2 * globRNG::uniform(globRNG::gen) - 1) * dmax;
      //std::cout << sx << " " << sy<< " " << sz << std::endl;

      vlist.move(particles, p, sx, sy, sz);
      //particles.applyBoundCond(p);

      // Compute epot_new - epot_old.
      delta_epot += potential.compute_Ui(particles, p, vlist);

      //draw a random number in (0,1]
      double rcomp = globRNG::uniform(globRNG::gen);

      if (rcomp > exp(-beta * delta_epot) and delta_epot > 0) //rejection
      {
        vlist.move(particles, p, -sx, -sy, -sz);
      }
    }//end loop over particles
  }//end loop over n_sweeps
}


void BoxSingleDisplace(Particles &BoxA, Particles &BoxB, Potential::WCA_U &potential, double dmax, double temperature, VerletList3D &vlistA, VerletList3D &vlistB, double &d_accRate)
{
  double boxRand = globRNG::uniform(globRNG::gen);

  if(boxRand < 0.5) // Displace a random particle in Box A.
  {
    SingleDisplace(BoxA, potential, dmax, temperature, vlistA, d_accRate);
  }

  else // Displace a random particle in Box B.
  {
    SingleDisplace(BoxB, potential, dmax, temperature, vlistB, d_accRate);
  }
}

void SingleDisplace(Particles &particles, Potential::WCA_U &potential, double dmax, double temperature, VerletList3D &vlist, double &d_accRate)
{
  //assume kb = 1;
  double beta = 1.0 / temperature;

  //choose a random particle
  int p = static_cast<int>(globRNG::uniform(globRNG::gen)*particles.n_particles_);

  //compute the old energy
  double delta_epot = -potential.compute_Ui(particles, p,vlist);

  //shift the particle position
  double sx = (2 * globRNG::uniform(globRNG::gen) - 1) * dmax;
  double sy = (2 * globRNG::uniform(globRNG::gen) - 1) * dmax;
  double sz = (2 * globRNG::uniform(globRNG::gen) - 1) * dmax;

  vlist.move(particles, p, sx, sy, sz);
  d_accRate++;

  //compute epot_trial-epot_new
  delta_epot += potential.compute_Ui(particles, p, vlist);

  //draw a random number in (0,1]
  double rcomp = globRNG::uniform(globRNG::gen);

  if (rcomp > exp(-beta * delta_epot) and delta_epot > 0) //rejection
  {
    vlist.move(particles, p, -sx, -sy, -sz);
    d_accRate--;
  }
}

};
