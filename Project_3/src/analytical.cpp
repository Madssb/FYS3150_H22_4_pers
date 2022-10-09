// Contains definitin of analytical solution
#include "analytical.hpp"
#include <cmath>
#include <armadillo>

arma::vec analytical_solution(arma::vec r_in, arma::vec v_in, double q,
                              double B_ext, double m, double V_ext,
                              double d, double t)
{
  arma::vec r = arma::vec(3);

  double x_0 = r_in(0);
  double z_0 = r_in(2);
  double v_0 = v_in(1);

  double w_0 = q * B_ext / m;
  double w_z_sqrd = 2. * q * V_ext / (m * std::pow(d, 2.0));
  double w_z = std::sqrt(w_z_sqrd);

  double w_plus = (w_0 + std::sqrt(std::pow(w_0, 2.0) - 2. * w_z_sqrd)) / 2.;
  double w_min = (w_0 - std::sqrt(std::pow(w_0, 2.0) - 2. * w_z_sqrd)) / 2.;

  double A_plus = (v_0 + w_min * x_0) / (w_min - w_plus);
  double A_min = - (v_0 + w_plus * x_0) / (w_min - w_plus);

  r(0) = A_plus * std::cos(w_plus * t) + A_min * std::cos(w_min * t);
  r(1) = - (A_plus * std::sin(w_plus * t) + A_min * std::sin(w_min * t));
  r(2) = z_0 * std::cos(w_z * t);

  return r;
}
