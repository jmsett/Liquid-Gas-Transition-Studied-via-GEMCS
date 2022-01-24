#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <limits>
#include "global_random.h"

class Particles
{

public:

  // Number particles, and Box-volume.
  int n_particles_;
  double BoxV_;

  // Box dimensions. (We will work with a cubic Box)
  double l_;
  double lx_;
  double ly_;
  double lz_;

  // Define space-coordinate vectors (3D System)
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<double> z_;

  Particles(int n, double Volume)
  {
    BoxV_ = Volume;
    n_particles_ = n;

    l_ = pow(BoxV_,1.0/3.0);
    lx_ = l_;
    ly_ = l_;
    lz_ = l_;

    x_.resize(n_particles_, 0.0);
    y_.resize(n_particles_, 0.0);
    z_.resize(n_particles_, 0.0);
  }

  Particles(std::string filename)
  {
    std::ifstream infile(filename.c_str());

    // checks if the file if open successfully
    if (!infile)
    {
        std::cout << "error: file not found" << std::endl;
        return;
    }

    std::string line;
    double x,y,z;
    int n;

    std::getline(infile, line);
    std::stringstream file(line);
    file >> BoxV_;
    BoxV_*=1.0;
    l_ = pow(BoxV_,1.0/3.0);
    lx_ = l_;
    ly_ = l_;
    lz_ = l_;

    while (std::getline(infile, line))// Read line by line
    {
        std::stringstream file(line);
        file >> n  >> x >> y >> z;
        x_.push_back(x);
        y_.push_back(y);
        z_.push_back(z);
    }
    n_particles_ = x_.size();
  }

  // Member functions (3D):

  void writeState(std::string filename)
  {
    std::ofstream outfile;
    outfile.open(filename.c_str());
    std::cout.precision(15);

    outfile << BoxV_ << std::endl;
    for (int i = 0; i < n_particles_ ; ++i)
    {
      outfile << i << ' ' << x_[i] << ' ' << y_[i] << ' ' << z_[i] << std::endl;
    }
    outfile.close();
  }

  void removePart(int i)
  {
    x_[i] = x_[n_particles_ - 1];
    x_.pop_back();
    y_[i] = y_[n_particles_ - 1];
    y_.pop_back();
    z_[i] = z_[n_particles_ - 1];
    z_.pop_back();

    n_particles_--;
  }

  void addParticle(double x, double y, double z)
  {
    x_.push_back(x);
    y_.push_back(y);
    z_.push_back(z);

    n_particles_++;
  }

  void setRandomPositions()
  {
    for (int i = 0; i < n_particles_; ++i)
    {
      x_[i] = globRNG::uniform(globRNG::gen) * lx_;
      y_[i] = globRNG::uniform(globRNG::gen) * ly_;
      z_[i] = globRNG::uniform(globRNG::gen) * lz_;
    }
  }

  void minimalImage(double &dx, double &dy, double &dz)
  {
    if (dx > 0.5 *  lx_)
    {
      dx -= lx_;
    }
    if (dx < 0.5 * -lx_)
    {
      dx += lx_;
    }

    if (dy > 0.5 *  ly_)
    {
      dy -= ly_;
    }
    if (dy < 0.5 * -ly_)
    {
      dy += ly_;
    }

    if (dz > 0.5 *  lz_)
    {
      dz -= lz_;
    }
    if (dz < 0.5 * -lz_)
    {
      dz += lz_;
    }
  }

  void applyBoundCond()
  {
    for(int i = 0; i < n_particles_; ++i)
    {
      if(x_[i] <= 0)  x_[i] += lx_;
      if(x_[i] > lx_) x_[i] -= lx_;

      if(y_[i] <= 0)  y_[i] += ly_;
      if(y_[i] > ly_) y_[i] -= ly_;

      if(z_[i] <= 0)  z_[i] += lz_;
      if(z_[i] > lz_) z_[i] -= lz_;
    }
  }

  void applyBoundCond(int i)
  {
    if(x_[i] <= 0)  x_[i] += lx_;
    if(x_[i] > lx_) x_[i] -= lx_;

    if(y_[i] <= 0)  y_[i] += ly_;
    if(y_[i] > ly_) y_[i] -= ly_;

    if(z_[i] <= 0)  z_[i] += lz_;
    if(z_[i] > lz_) z_[i] -= lz_;
  }
};
