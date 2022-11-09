/*
Source file for the functions to calculate expected values
*/

#include "expectation_values.hpp"
#include <map>
#include <iostream>
#include <cmath>

// Calculate the expected energy
void expected_energy(std::vector<double>& expEpsilon,
                     std::vector<double>& expEpsilon_squared,
                     std::vector<double> epsilonEnergies,
                     double initEnergy, double T)
{
  double N = epsilonEnergies.size();

  double expectedEnergy = initEnergy;
  double expectedEnergy_squared = initEnergy * initEnergy;

  expEpsilon.push_back(initEnergy);
  expEpsilon_squared.push_back(initEnergy * initEnergy);

  for (int n = 1; n <= N; n++)
  {
    int E = epsilonEnergies[n];
    expectedEnergy += E;
    expectedEnergy_squared += (E * E);

    expEpsilon.push_back(expectedEnergy / n);
    expEpsilon_squared.push_back(expectedEnergy_squared / n);
  }
}
