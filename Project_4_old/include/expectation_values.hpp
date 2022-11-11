/*
Header file for the functions to calculate the expected values
*/

#ifndef __Expectation_Values_hpp__
#define __Expectation_Values_hpp__
#include <vector>

void expected_values(std::vector<double>& expEpsilon,
                     std::vector<double>& expEpsilon_squared,
                     std::vector<double> epsilonEnergies,
                     std::vector<double> magnPerSpin,
                     std::vector<double>& expMagn,
                     std::vector<double>& expMagn_squared,
                     double initMagn, double initEnergy);

#endif
