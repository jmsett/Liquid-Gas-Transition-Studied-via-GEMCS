#pragma once

#include <iostream>
#include "celllist.h"

class VerletList3D
{
public:

    int nparticles_;
    double cutoff_;
    double skin_;
    double skin_sq_;
    double cutoff_tot_;
    double cutoff_tot_sq_;
    double container_l_;
    double lcell_;

    int nverlet_;

    std::vector<int> head_;
    std::vector<int> list_;

    std::vector<double> dx_;
    std::vector<double> dy_;
    std::vector<double> dz_;

private:
    CellList * celllist_;

public:
    VerletList3D(int n_particles, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, double container_l, double cutoff)
    {
        nparticles_ = n_particles;
        cutoff_ = cutoff;
        skin_ = cutoff * 0.3;

        skin_sq_ = skin_ * skin_;
        cutoff_tot_ = cutoff_ + skin_;
        cutoff_tot_sq_ = cutoff_tot_ * cutoff_tot_;
        container_l_ = container_l;

        /*
        double container_v = std::pow(container_l, 3);
        int nneighapprox (static_cast<int> ((n_particles/container_v)*(4./3.) * cutoff_tot_sq_ * cutoff_tot_ * 3.1416));
        if (nneighapprox < 2) nneighapprox = 4;
        nverlet_ = n_particles * nneighapprox * 2;
        */
        // Worst case: Starting from small density, but high density clusters form.
        nverlet_ = n_particles * n_particles;

        head_.resize(n_particles + 1, 0);
        list_.resize(nverlet_, 0);
        dx_.resize(n_particles, 0.0);
        dy_.resize(n_particles, 0.0);
        dz_.resize(n_particles, 0.0);
        celllist_ = new CellList(n_particles, x, y, z, container_l, cutoff_, skin_);
        lcell_ = celllist_->lcell_;
        update(x, y, z);
    }

    void reset()
    {
      // Reset vector values
      head_.clear();
      list_.clear();
      dx_.clear();
      dy_.clear();
      dz_.clear();
      //std::fill(head_.begin(), head_.end(), 0);
      //std::fill(list_.begin(), list_.end(), 0);
      //std::fill(dx_.begin(), dx_.end(), 0.0);
      //std::fill(dy_.begin(), dy_.end(), 0.0);
      //std::fill(dz_.begin(), dz_.end(), 0.0);

      // reset member variables
      nparticles_=0;
      cutoff_=0.0;
      skin_=0.0;
      skin_sq_=0.0;
      cutoff_tot_=0.0;
      cutoff_tot_sq_=0.0;
      container_l_=0.0;
      lcell_=0.0;
      //nverlet_=0;

      //reset CellList
      celllist_->reset();
    }

    void redefine(int n_particles, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, double container_l, double cutoff)
    {
        nparticles_ = n_particles;
        cutoff_ = cutoff;
        skin_ = cutoff * 0.3;

        skin_sq_ = skin_ * skin_;
        cutoff_tot_ = cutoff_ + skin_;
        cutoff_tot_sq_ = cutoff_tot_ * cutoff_tot_;
        container_l_ = container_l;

        /*
        double container_v = std::pow(container_l, 3);
        int nneighapprox (static_cast<int> ((n_particles/container_v)*(4./3.) * cutoff_tot_sq_ * cutoff_tot_ * 3.1416));
        if (nneighapprox < 2) nneighapprox = 4;
        nverlet_ = n_particles * nneighapprox * 2;
        */
        nverlet_ = n_particles * n_particles;

        head_.resize(n_particles + 1, 0);
        list_.resize(nverlet_, 0);
        dx_.resize(n_particles, 0.0);
        dy_.resize(n_particles, 0.0);
        dz_.resize(n_particles, 0.0);

        celllist_->redefine(n_particles, x, y, z, container_l, cutoff_, skin_);
        lcell_ = celllist_->lcell_;
        update(x, y, z);
    }

    void v_update(double container_l, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)
    {
        container_l_ = container_l;

        double container_v = std::pow(container_l, 3);
        int nneighapprox (static_cast<int> ((nparticles_/container_v)*(4./3.) * cutoff_tot_sq_ * cutoff_tot_ * 3.1416));
        if (nneighapprox < 2) nneighapprox = 4;
        nverlet_ = nparticles_ * nneighapprox * 2;

        list_.resize(nverlet_, 0);

        delete celllist_;
        celllist_ = new CellList(nparticles_, x, y, z, container_l, cutoff_, skin_);
        update(x, y, z);
    }

    void n_update(int n_particles, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)
    {
      nparticles_ = n_particles;

      double container_v = std::pow(container_l_, 3);
      int nneighapprox (static_cast<int> ((n_particles/container_v)*(4./3.) * cutoff_tot_sq_ * cutoff_tot_ * 3.1416));
      if (nneighapprox < 2) nneighapprox = 4;
      nverlet_ = n_particles * nneighapprox * 2;

      head_.resize(n_particles + 1, 0);
      list_.resize(nverlet_, 0);
      dx_.resize(n_particles, 0.0);
      dy_.resize(n_particles, 0.0);
      dz_.resize(n_particles, 0.0);

      //celllist_->n_update(n_particles);
      delete celllist_;
      celllist_ = new CellList(nparticles_, x, y, z, container_l_, cutoff_, skin_);
      update(x, y, z);
    }

    void update(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)
    {
        celllist_->update(nparticles_, x, y, z);
        int pcx, pcy, pcz; // Particle Cell
        int ccx, ccy, ccz; // Current Cell
        int sx, sy, sz;    // Stencil index

        int nx = celllist_->nx_;

        int ravel_xy = nx*nx;

        double x1, y1, z1, x2, y2, z2, dx12, dy12, dz12;

        int p2 = 0;
        int verlet_list_offset = 0;

        //reseting the particle traveled distance bookeeping

        std::fill(dx_.begin(), dx_.end(), 0.0);
        std::fill(dy_.begin(), dy_.end(), 0.0);
        std::fill(dz_.begin(), dz_.end(), 0.0);

        for (int p1 = 0; p1 < nparticles_; ++p1)
        {
            x1 = x[p1];
            y1 = y[p1];
            z1 = z[p1];

            pcx = static_cast<int>(x[p1]/lcell_);
            pcy = static_cast<int>(y[p1]/lcell_);
            pcz = static_cast<int>(z[p1]/lcell_);

            for (int s = 0; s < celllist_->nstencil_; ++s) // iterate over all the cells pointed by the stencil
            {
                sx = celllist_->stencil_x_[s];
                sy = celllist_->stencil_y_[s];
                sz = celllist_->stencil_z_[s];

                ccx = sx + pcx;
                ccy = sy + pcy;
                ccz = sz + pcz;

                //-------------------------
                if (ccx < 0)   ccx += nx;
                if (ccx >= nx) ccx -= nx;
                //-------------------------
                if (ccy < 0)   ccy += nx;
                if (ccy >= nx) ccy -= nx;
                //--------------------------
                if (ccz < 0)   ccz += nx;
                if (ccz >= nx) ccz -= nx;
                //--------------------------

                p2 = celllist_->head_cells_[ccx + ccy * nx + ccz * ravel_xy];

                while (p2 >= 0)
                {
                    x2 = x[p2];
                    y2 = y[p2];
                    z2 = z[p2];

                    dx12 = x1 - x2;
                    dy12 = y1 - y2;
                    dz12 = z1 - z2;

                    //------Minimal Image Convention
                    if (dx12 <= -0.5 * container_l_) dx12 += container_l_;
                    if (dx12 >   0.5 * container_l_) dx12 -= container_l_;
                    if (dy12 <= -0.5 * container_l_) dy12 += container_l_;
                    if (dy12 >   0.5 * container_l_) dy12 -= container_l_;
                    if (dz12 <= -0.5 * container_l_) dz12 += container_l_;
                    if (dz12 >   0.5 * container_l_) dz12 -= container_l_;
                    //------------------------------

                    if (dx12 * dx12 + dy12 * dy12 + dz12 * dz12 <= cutoff_tot_sq_ and p1 != p2)
                    {
                        list_[verlet_list_offset] = p2;
                        verlet_list_offset ++;

                        //comment out range checking for speedup
                        if(verlet_list_offset >= nverlet_)
                        {
                          std::cout << std::endl;
                          std::cout << "verlet_list overflow" << std::endl;
                          //std::cout << "(p1,p2): (" << p1 << "," << p2 << ") nverlet: " << nverlet_ << " list_size: "<<list_.size() << " offset: " << verlet_list_offset << " list[offset]: " << list_[verlet_list_offset] << std::endl;
                          //exit(0);
                        }

                    }
                    p2 = celllist_->list_cells_[p2];
                }
            }
            head_[p1 + 1] = verlet_list_offset;
        }
    }

    inline void move (Particles &particles, int n, double dx, double dy, double dz)
{
    dx_[n] += dx;
    dy_[n] += dy;
    dz_[n] += dz;

    particles.x_[n] += dx;
    particles.y_[n] += dy;
    particles.z_[n] += dz;

    particles.applyBoundCond(n);

    if (dx_[n] * dx_[n] + dy_[n] * dy_[n] + dz_[n] * dz_[n] >= 0.25 * skin_sq_)
    {
        update(particles.x_, particles.y_, particles.z_);
    }
}

~VerletList3D()
{
    delete celllist_;
}

};
