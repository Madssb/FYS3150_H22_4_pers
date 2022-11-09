/*
Header file for the functions to calculate the expected values
*/

#ifndef __Expectation_Values_hpp__
#define __Expectation_Values_hpp__
#include <vector>

void expected_energy(std::vector<double>& expEpsilon,
                     std::vector<double>& expEpsilon_squared,
                     std::vector<double> epsilonEnergies,
                     double initEnergy, double T);

std::vector<double> expected_magnetization(std::vector<double> magnetization,
                                           double T);

#endif
