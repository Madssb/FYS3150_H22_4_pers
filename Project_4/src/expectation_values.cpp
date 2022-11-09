/*
Source file for the functions to calculate expected values
*/

#include "expectation_values.hpp"
#include <map>
#include <iostream>
#include <cmath>

// Calculate the expected energy
std::vector<double> expected_energy(std::vector<double> epsilonEnergies,
                                    double initEnergy, double T)
{
  double N = epsilonEnergies.size();
  std::vector<double> expEpsilon;
  double expectedEnergy = initEnergy;
  expEpsilon.push_back(initEnergy);

  // Partition funct. for 2x2 case
  double Z = 2 * (std::cosh(8 / T) + 3);

  for (int n = 1; n <= N; n++)
  {
    int E = epsilonEnergies[n];
    expectedEnergy += E;
    expEpsilon.push_back(expectedEnergy / n);
  }

  return expEpsilon;
}
