#pragma once
#include <random>
//another method would be to implement the singleton sheme...
namespace globRNG
{
  std::random_device rd;
  std::mt19937 gen(rd());

  //declare global random number distributions
  std::uniform_real_distribution<double> uniform(0, 1);
}
