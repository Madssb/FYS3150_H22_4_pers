/*
This is the main code for solving the time-dependent Schrödinger equation
*/
#include "matrix.hpp"

using namespace std;
using namespace arma;

int main(int argc, const char* argv[])
{
  const string configuration_in = argv[2];
  const string potential_in = argv[3];

  string config_file = configuration_in + ".txt";
  string potential = potential_in + ".txt";

  arma::rowvec config;
  config.load(config_file, arma::raw_ascii);

  double h_in = config(0);
  double dt_in = config(1);
  double T_in = config(2);
  double x_c = config(3);
  double sigma_x = config(4);
  double p_x = config(5);
  double y_c = config(6);
  double sigma_y = config(7);
  double p_y = config(8);
  double V_0 = config(9);



  Matrix test = Matrix(h_in, dt_in, T_in);
  test.fill_matrices();
  test.set_initial_state(x_c, sigma_x, p_x, y_c, sigma_y, p_y);
  test.set_potential(potential, V_0);
  // test.solve();

  return 0;
}
