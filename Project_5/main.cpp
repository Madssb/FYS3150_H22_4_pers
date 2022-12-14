/*
This is the main code for solving the time-dependent Schrödinger equation
*/
#include "matrix.hpp"

using namespace std;
using namespace arma;

int main(int argc, const char* argv[])
{
  const string configuration = argv[2];
  const string potential = argv[3];
  const string filename = argv[4];

  string out_filename = filename + ".bin";

  arma::rowvec config;
  config.load(configuration, arma::raw_ascii);

  double h_in = config(0);
  double dt_in = config(1);
  double T_in = config(2);
  double x_c = config(3);
  double sigma_x = config(4);
  double p_x = config(5);
  double y_c = config(6);
  double sigma_y = config(7);
  double p_y = config(8);



  Matrix matrix = Matrix(h_in, dt_in, T_in);
  matrix.set_potential(potential);
  matrix.fill_matrices();
  matrix.set_initial_state(x_c, sigma_x, p_x, y_c, sigma_y, p_y);
  matrix.solve();
  matrix.S.save(out_filename);

  return 0;
}
