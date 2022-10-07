// Contains definitin of analytical solution
#include "analytical_z.hpp"
#include <cmath>

double analytical_solution_z(double q, double V_0, double m, double d, double t)
{
  double w_z_squared = 2 * q * V_0 / (m * pow(d, 2.));
  double w_z = sqrt(w_z_squared);
  double z_0 = 1.;
  double z = z_0 * cos(w_z * t);

  return z;
}
