#pragma once

#include <iostream>
#include <fstream>

#include "particles.h"
#include "montecarlo.h"
#include "progresBar.h"

namespace CorrelationsFunctions
{
  std::vector<double> gr(Particles&, double);

  // nbins is th number of "histogram"
  std::vector<double> gr(Particles &particles, double nbins)
  {
    //choose the smallest dimension of the container
    double length = particles.lx_;
    if( length > particles.ly_) length = particles.ly_;

    //half of container lenght will the considered max distance
    length /= 2.0;
    double bin_length = length/(nbins);

    //number_density
    double number_density = (particles.n_particles_)/(particles.lx_*particles.ly_*particles.lz_);

    std::vector<double> bins(nbins, 0.0);
    int current_bin;
    double dx, dy, dz;
    double dr;

    for (int i = 0; i < particles.n_particles_; ++i)
    {
      for (int j = 0; j < i ; ++j)
      {
        dx = particles.x_[i]-particles.x_[j];
        dy = particles.y_[i]-particles.y_[j];
        dz = particles.z_[i]-particles.z_[j];

        particles.minimalImage(dx,dy,dz);

        dr = sqrt(dx * dx + dy * dy + dz * dz);
        current_bin = floor(nbins * dr/length);
        if(current_bin < nbins)
        {
          bins[current_bin] = bins[current_bin] + 2.0;
        }
      }
    }

    //nomrmalization and output
    double shell_volume;
    for (int i = 0; i < nbins; ++i)
    {
      shell_volume =  (4.0/3.0)*M_PI * (pow((i+1) * bin_length, 3) - pow(i * bin_length, 3));
      bins[i] /=  shell_volume * particles.n_particles_ * number_density;

    }
    //bins[0] = 0;
    return bins;
  }

  void calculate_gr(int nbin, int n_avgs, Particles particles, Potential::WCA_U &potential, double dmax, double temperature, VerletList3D& vlist)
  {
    progresBar bar(n_avgs, "Correlation");

    std::vector<double> gr_m;
    std::vector<double> gr_m_final(nbin, 0.0);

    for (int i = 0; i < n_avgs; i++)
    {
      bar.update();
      gr_m = gr(particles, nbin);
      Montecarlo::Canonical2(particles, potential, 10, dmax, temperature, vlist);

      for (int i = 0; i < nbin; i++)
      {
         gr_m_final[i] += gr_m[i];
      }
    }
    bar.end();

    //average and write gr to file
    std::ofstream gr_file("gr.dat");
    for (int i = 0; i < nbin; i++)
    {
      gr_file << i * 0.5 * particles.lx_/nbin << ' ' << gr_m_final[i]/n_avgs <<  std::endl;
    }
  }

}
