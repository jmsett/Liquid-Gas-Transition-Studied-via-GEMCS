#pragma once

#include <cmath>
#include <vector>

class CellList
{
public:

    long ncells_;
    long nx_;
    double lcell_;
    int n_cell_dim;

    int sh_;
    int sl_;

    double volume_;

    int ravel_xy_;
    int nstencil_;

    std::vector<int> head_cells_;
    std::vector<int> list_cells_;

    std::vector<int> stencil_x_;
    std::vector<int> stencil_y_;
    std::vector<int> stencil_z_;

    void reset()
    {
      stencil_x_.clear();
      stencil_y_.clear();
      stencil_z_.clear();
      head_cells_.clear();
      list_cells_.clear();

      ncells_=0;
      nx_=0;
      lcell_=0.0;
      n_cell_dim=0;

      sh_=0;
      sl_=0;

      volume_=0.0;

      ravel_xy_=0;
      nstencil_=0;
    }

    //constructor for 3D only for cubic boxes !!!!!!!
    CellList(int n_particles, std::vector<double> &particles_x, std::vector<double> &particles_y, std::vector<double> &particles_z, double container_l, double cutoff, double skin)
    {
        lcell_ = cutoff * 0.5;
        nx_ = std::ceil(container_l / lcell_);
        lcell_ = container_l/nx_;

        ncells_ = nx_ * nx_ * nx_;
        ravel_xy_ = nx_ * nx_;

        head_cells_.resize(ncells_, 0.0);
        list_cells_.resize(n_particles, 0.0);

        nstencil_ = 0;

        double cutoff_tot = cutoff + skin;
        double cutoff_tot_sq = cutoff_tot * cutoff_tot;

        sh_ = ceil ((cutoff_tot) / (lcell_));
        sl_ = -sh_;

        if (2 * sh_ + 1 >= nx_)
        {
            sh_ = nx_ / 2;
            sl_ =-nx_ / 2 + (nx_ + 1) % 2;
        }

        for (int i = sl_; i <= sh_; ++i) {
            for ( int j = sl_; j <= sh_; ++j) {
                for ( int k = sl_; k <= sh_; ++k) {
                    if (cellDistanceSquared(i,j,k) <= cutoff_tot_sq) {
                        stencil_x_.push_back(i);
                        stencil_y_.push_back(j);
                        stencil_z_.push_back(k);

                    }
                }
            }
        }
        nstencil_ = static_cast<int>(stencil_x_.size());

        int bin_nr = 0;
        int i, j, k;

        for (int i = 0; i < ncells_; ++i)
        {
            head_cells_[i] = -1;
        }

        stencil_x_.shrink_to_fit();
        stencil_y_.shrink_to_fit();
        stencil_z_.shrink_to_fit();

        for (int n = 0; n < n_particles; ++n)
        {
            i = static_cast<int>(particles_x[n]/lcell_);
            j = static_cast<int>(particles_y[n]/lcell_);
            k = static_cast<int>(particles_z[n]/lcell_);

            bin_nr = i + j * nx_ + k * ravel_xy_;

            list_cells_[n] = head_cells_[bin_nr];
            head_cells_[bin_nr] = n;
        }
        return;
    }

    void redefine(int n_particles, std::vector<double> &particles_x, std::vector<double> &particles_y, std::vector<double> &particles_z, double container_l, double cutoff, double skin)
    {
        lcell_ = cutoff * 0.5;
        nx_ = std::ceil(container_l / lcell_);
        lcell_ = container_l/nx_;

        ncells_ = nx_ * nx_ * nx_;
        ravel_xy_ = nx_ * nx_;

        head_cells_.resize(ncells_, 0.0);
        list_cells_.resize(n_particles, 0.0);

        nstencil_ = 0;

        double cutoff_tot = cutoff + skin;
        double cutoff_tot_sq = cutoff_tot * cutoff_tot;

        sh_ = ceil ((cutoff_tot) / (lcell_));
        sl_ = -sh_;

        if (2 * sh_ + 1 >= nx_)
        {
            sh_ = nx_ / 2;
            sl_ =-nx_ / 2 + (nx_ + 1) % 2;
        }

        for (int i = sl_; i <= sh_; ++i) {
            for ( int j = sl_; j <= sh_; ++j) {
                for ( int k = sl_; k <= sh_; ++k) {
                    if (cellDistanceSquared(i,j,k) <= cutoff_tot_sq) {
                        stencil_x_.push_back(i);
                        stencil_y_.push_back(j);
                        stencil_z_.push_back(k);

                    }
                }
            }
        }
        nstencil_ = static_cast<int>(stencil_x_.size());

        int bin_nr = 0;
        int i, j, k;

        for (int i = 0; i < ncells_; ++i)
        {
            head_cells_[i] = -1;
        }

        stencil_x_.shrink_to_fit();
        stencil_y_.shrink_to_fit();
        stencil_z_.shrink_to_fit();

        for (int n = 0; n < n_particles; ++n)
        {
            i = static_cast<int>(particles_x[n]/lcell_);
            j = static_cast<int>(particles_y[n]/lcell_);
            k = static_cast<int>(particles_z[n]/lcell_);

            bin_nr = i + j * nx_ + k * ravel_xy_;

            list_cells_[n] = head_cells_[bin_nr];
            head_cells_[bin_nr] = n;
        }
        return;
    }

    //constructor for 2D
    CellList(int n_particles, std::vector<double> &particles_x, std::vector<double> &particles_y, double container_l, double cutoff, double skin)
    {
        lcell_ = cutoff * 0.5;
        nx_ = ceil(container_l / lcell_);
        lcell_ = container_l/nx_;

        ncells_ = nx_ * nx_;
        ravel_xy_ = nx_;

        head_cells_.resize(ncells_, 0.0);
        list_cells_.resize(n_particles, 0.0);

        double cutoff_tot = cutoff + skin;
        double cutoff_tot_sq = cutoff_tot * cutoff_tot;

        sh_ = ceil (cutoff_tot / lcell_);
        sl_ = -sh_;

        if (2 * sh_ + 1 >= nx_)
        {
            sh_ = nx_ / 2;
            sl_ =-nx_ / 2 + (nx_ + 1) % 2;
        }

        for (int i = sl_; i <= sh_; i++) {
            for ( int j = sl_; j <= sh_; j++) {

                if (cellDistanceSquared(i,j) <= cutoff_tot_sq) {
                    stencil_x_.push_back(i);
                    stencil_y_.push_back(j);
                }
            }
        }

        nstencil_ = static_cast<int>(stencil_x_.size());

        stencil_x_.shrink_to_fit();
        stencil_y_.shrink_to_fit();

        int bin_nr = 0;
        int i, j;

        for (int i = 0; i < ncells_; ++i)
        {
            head_cells_[i] = -1;
        }

        for (int n = 0; n < n_particles; ++n)
        {
            i = std::floor(particles_x[n]/lcell_);
            j = std::floor(particles_y[n]/lcell_);

            bin_nr = i + j * nx_;

            list_cells_[n] = head_cells_[bin_nr];
            head_cells_[bin_nr] = n;
        }
    }

    // -------- Added Code----------------
    void n_update(int n_particles)
    {
      list_cells_.resize(n_particles, 0.0);
    }
    // -----------------------------------

    void update(int n_particles, std::vector<double> &particles_x, std::vector<double> &particles_y, std::vector<double> &particles_z)
    {
        int bin_nr = 0;
        int i, j, k;

        for (int n = 0; n < ncells_; ++n)
        {
            head_cells_[n] = -1;
        }

        for (int n = 0; n < n_particles; ++n)
        {
            i = static_cast<int>(particles_x[n]/lcell_);
            j = static_cast<int>(particles_y[n]/lcell_);
            k = static_cast<int>(particles_z[n]/lcell_);

            bin_nr = i + j * nx_ + k * ravel_xy_;
            list_cells_[n] = head_cells_[bin_nr];
            head_cells_[bin_nr] = n;
        }
    }

    void update(int n_particles,std::vector<double> &particles_x, std::vector<double> &particles_y)
    {
        int bin_nr = 0;
        int i, j;

        for (int n = 0; n < ncells_; ++n)
        {
            head_cells_[n] = -1;
        }

        for (int n = 0; n < n_particles; ++n)
        {
            i = std::floor(particles_x[n]/lcell_);
            j = std::floor(particles_y[n]/lcell_);

            bin_nr = i + j * nx_;
            list_cells_[n] = head_cells_[bin_nr];
            head_cells_[bin_nr] = n;
        }
    }

private:

    /* ----------------------------------------------------------------------
     *   compute closest distance between central bin (0,0,0) and bin (i,j,k)
     * ------------------------------------------------------------------------- */

    double cellDistanceSquared(int i, int j, int k)
    {
        double dx, dy,dz;
        //-------------------------------x
        if( i < 0)      dx = this->lcell_*(i+1);
        else if(i == 0) dx = 0.0;
        else            dx = this->lcell_*(i-1);
        //-------------------------------y
        if ( j < 0)     dy = this->lcell_*(j+1);
        else if(j == 0) dy = 0.0;
        else            dy = this->lcell_*(j-1);
        //-------------------------------z
        if ( k < 0)     dz = this->lcell_*(k+1);
        else if(k == 0) dz = 0.0;
        else            dz = this->lcell_*(k-1);
        return dx*dx+dy*dy+dz*dz;
    }

    double cellDistanceSquared(int i, int j)
    {
        double dx, dy;
        //-------------------------------x
        if( i < 0)      dx = this->lcell_*(i+1);
        else if(i == 0) dx = 0.0;
        else            dx = this->lcell_*(i-1);
        //-------------------------------y
        if ( j < 0)     dy = this->lcell_*(j+1);
        else if(j == 0) dy = 0.0;
        else            dy = this->lcell_*(j-1);

        return dx*dx+dy*dy;
    }
};
