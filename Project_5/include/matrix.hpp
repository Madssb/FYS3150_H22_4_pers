/*
This file contains the definition of the class Matrix
*/

#include <armadillo>

class Matrix
{
  public:

    // Member variables
    double h, dt, T;        // Spatial and time step length and total time
    int M, N;               // No. of points along x, y and time axes

    arma::cx_double r;      // Predefined constant, r = i∆t/h^2

    arma::sp_cx_mat A, B;   // Crank-Nicolson matrices
    arma::cx_mat U;         // State matrix
    arma::mat V;            // Potential matrix

    arma::cx_cube S;        // Storing states

    // Constructor
    Matrix(double h_in, double dt_in, double T_in);

    // Takes indices (i, j) and returns the corresponding k index
    int index_k(int i, int j);

    // Set potential from file
    void set_potential(std::string filename, double V_0);

    // Fill matrices A and B
    void fill_matrices();

};
