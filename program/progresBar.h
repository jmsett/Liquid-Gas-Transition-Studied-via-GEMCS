#pragma once
#include <iostream>

class progresBar
{
  double delta_;
  double state_;
  std::string process_;
  std::string subprocess_;

public:

  progresBar(int steps, std::string process)
  {
    delta_ = 100.0 / ((steps)*1.0);
    state_ = 0.0;
    process_ = process;
  }

  progresBar(int steps, std::string process, std::string subprocess)
  {
    delta_ = 100.0 / ((steps)*1.0);
    state_ = 0.0;
    process_ = process;
    subprocess_ = subprocess;
  }

  void name(std::string newProcess)
  {
    process_ = newProcess;
  }

  void subname(std::string newProcess)
  {
    subprocess_ = newProcess;
  }

  void update()
  {
    state_ += delta_;
    std::cout << " " << process_ <<" [" << std::floor(state_) << "%] "
    << "                                       \r" << std::flush;
  }

  void update(double intensV)
  {
    state_ += delta_;
    std::cout << " " << process_ <<" [" << std::floor(state_) << "%] "
    << intensV << "                                       \r" << std::flush;
  }

  void update(int c1, int c2, int c3)
  {
    state_ += delta_;
    std::cout << " " << process_ << " [" << std::floor(state_) << "%]"
    << "  :  " << subprocess_ << " (" << c1 << ", " << c2 << ", " << c3 << ")"
    << "                      \r" << std::flush;
  }

  void end()
  {
    std::cout << "                                                                               \r" << std::flush;
    std::cout << process_ << " ... Done" << std::endl;
  }

};
