/*
Source file for the functions to calculate expected values
*/

#include "expectation_values.hpp"
#include <map>
#include <iostream>
#include <cmath>

// Calculate the expected energy
void expected_values(std::vector<double>& expEpsilon,
                     std::vector<double>& expEpsilon_squared,
                     std::vector<double> epsilonEnergies,
                     std::vector<double> magnPerSpin,
                     std::vector<double>& expMagn,
                     std::vector<double>& expMagn_squared,
                     double initMagn, double initEnergy)
{
  double N = epsilonEnergies.size();

  double expectedEnergy = initEnergy;
  double expectedEnergy_squared = initEnergy * initEnergy;

  double expectedMagn = initMagn;
  double expectedMagn_squared = initMagn * initMagn;

  expEpsilon.push_back(initEnergy);
  expEpsilon_squared.push_back(initEnergy * initEnergy);
  expMagn.push_back(abs(initMagn));
  expMagn_squared.push_back(initMagn * initMagn);

  for (int n = 1; n <= N; n++)
  {
    double E = epsilonEnergies[n];
    double M = magnPerSpin[n];

    expectedEnergy += E;
    expectedEnergy_squared += (E * E);

    expectedMagn += abs(M);
    expectedMagn_squared += (M * M);

    expEpsilon.push_back(expectedEnergy / n);
    expEpsilon_squared.push_back(expectedEnergy_squared / n);

    expMagn.push_back(expectedMagn / n);
    expMagn_squared.push_back(expectedMagn_squared / n);
  }
}
