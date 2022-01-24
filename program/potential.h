#pragma once

#include <cmath>
#include "particles.h"
#include "verletlist3D.h"

namespace Potential
{

class WCA_U
{
  public:

    //member variables
    double epsilon_;
    double cutoff_;
    double cutoff_sq_;
    double U_rc_;
    double dU_rc_;

    double utot6_;
    double utot12_;
    double U_rc_tot_;
    double dU_rc_tot_;
    double dU_r_tot_;


    // constructor
    WCA_U(double epsilon, double cutoff)
    {
      epsilon_ = epsilon;
      cutoff_ = cutoff;
      cutoff_sq_ = cutoff_ * cutoff_;

      double rc_1_inv = 1.0 / cutoff_;
      double rc_2_inv = 1.0 / cutoff_sq_;
      double rc_6_inv = rc_2_inv * rc_2_inv * rc_2_inv;

      U_rc_ = 4. * epsilon_ * rc_6_inv * (rc_6_inv - 1.0);
      dU_rc_ = -48. * epsilon_ * rc_1_inv * rc_6_inv * (rc_6_inv - 0.5);
    }

    inline void pot(double r)
    {
      double r_inv = 1.0 / r;
      double r_2_inv = r_inv * r_inv;
      double r_6_inv = r_2_inv * r_2_inv * r_2_inv;

      utot12_ += 4. * epsilon_ * r_6_inv * r_6_inv;
      utot6_ += 4. * epsilon_ * r_6_inv;
      U_rc_tot_ += U_rc_;
      dU_rc_tot_ += cutoff_ * dU_rc_;
      dU_r_tot_ += r * dU_rc_;
    }

    inline double pot_i(double r)
    {
      double r_inv = 1.0 / r;
      double r_2_inv = r_inv * r_inv;
      double r_6_inv = r_2_inv * r_2_inv * r_2_inv;

      return 4.*epsilon_*r_6_inv*(r_6_inv - 1.0) - U_rc_ - (r - cutoff_)*dU_rc_;
    }

    inline double force(double r)
    {
      double r_inv = 1.0 / r;
      double r_2_inv = r_inv * r_inv;
      double r_6_inv = r_2_inv * r_2_inv * r_2_inv;

      //  F = -d/dr (Uij)
      return (-1.0)*(-48.*epsilon_*r_6_inv*r_inv*(r_6_inv - 0.5) - dU_rc_);
    }

    double compute_P(double Temp, Particles &particles, VerletList3D &vlist){

      double sum = 0.0;

      double dx, dy, dz;
      double r_abs = 0.0;

      for (int i = 0; i < particles.n_particles_; ++i)
      {
        for (int j = vlist.head_[i]; j < vlist.head_[i+1]; ++j)
        {
          int p = vlist.list_[j];
          if (i < p) continue;

          dx = particles.x_[i] - particles.x_[p];
          dy = particles.y_[i] - particles.y_[p];
          dz = particles.z_[i] - particles.z_[p];

          particles.minimalImage(dx, dy, dz);

          r_abs = dx * dx + dy * dy + dz * dz;
          r_abs = sqrt(r_abs);

          if (r_abs < cutoff_)
          {
            sum += r_abs*force(r_abs);
          }
        }
      }
      return (particles.n_particles_*Temp + sum/3.0) / (1.0*particles.BoxV_);
    }


    double scaled_Utot(double l_old, double l_new)
    {
      double scaling = l_old/l_new;
      double scaling_inv = 1.0/scaling;
      double scaling_6 = scaling * scaling * scaling * scaling * scaling * scaling;
      double scaling_12 = scaling_6 * scaling_6;

      return (utot12_*scaling_12 - utot6_*scaling_6) - U_rc_tot_ - (dU_r_tot_*scaling_inv - dU_rc_tot_);
    }

    double compute_Utot(Particles &particles, VerletList3D &vlist)
    {
      utot12_ = 0.0;
      utot6_ = 0.0;
      U_rc_tot_ = 0.0;
      dU_rc_tot_ = 0.0;
      dU_r_tot_ = 0.0;

      double dx, dy, dz;
      double r_abs;

      for (int i = 0; i < particles.n_particles_; ++i)
      {
        for (int j = vlist.head_[i]; j < vlist.head_[i+1]; ++j)
        {
          int p = vlist.list_[j];
          if (i < p) continue;

          dx = particles.x_[i] - particles.x_[p];
          dy = particles.y_[i] - particles.y_[p];
          dz = particles.z_[i] - particles.z_[p];

          particles.minimalImage(dx, dy, dz);

          r_abs = dx * dx + dy * dy + dz * dz;
          r_abs = sqrt(r_abs);

          if (r_abs < cutoff_)
          {
            pot(r_abs);
          }
        }
      }
      return (utot12_ - utot6_) - U_rc_tot_ - (dU_r_tot_- dU_rc_tot_);
    }

    double compute_Ui(Particles &particles, int i, VerletList3D &vlist)
    {
      double u_tot = 0.0;
      double dx, dy, dz;
      double r_abs;

      for (int j = vlist.head_[i]; j < vlist.head_[i+1]; ++j)
      {
        int p = vlist.list_[j];

        dx = particles.x_[i] - particles.x_[p];
        dy = particles.y_[i] - particles.y_[p];
        dz = particles.z_[i] - particles.z_[p];

        particles.minimalImage(dx, dy, dz);

        r_abs = dx * dx + dy * dy + dz * dz;
        r_abs = sqrt(r_abs);

        if (r_abs < cutoff_)
        {
          u_tot += pot_i(r_abs);
        }
      }
      return u_tot;
    }

    double computeVirtual(Particles &particles, double x, double y, double z)
    {
      double u_tot = 0.0;
      double dx, dy, dz;
      double r_abs;

      for (int j = 0; j < particles.n_particles_; ++j)
      {
        dx = x - particles.x_[j];
        dy = y - particles.y_[j];
        dz = z - particles.z_[j];

        particles.minimalImage(dx, dy, dz);

        r_abs = dx * dx + dy * dy + dz * dz;
        r_abs = sqrt(r_abs);

        if (r_abs < cutoff_)
        {
          u_tot += pot_i(r_abs);
        }
      }
      return u_tot;
    }

}; // END Class
}  // END Namespace
