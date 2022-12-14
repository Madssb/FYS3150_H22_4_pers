/*
This file contains the definition of the class Matrix
*/

#include <armadillo>

class Matrix
{
  public:

    // Member variables
    double h, dt, T;                // Spatial and time step length and total time
    int M, N;                       // No. of points along x, y and time axes

    arma::cx_double r;              // Predefined constant, r = i∆t/h^2

    arma::sp_cx_mat A, B;           // Crank-Nicolson matrices
    arma::cx_mat U;                 // State matrix
    arma::cx_vec u_new, u_current;  // Column vectors for u^(n+1) and u^n
    arma::mat V;                    // Potential matrix

    arma::cx_cube S;                // Storing states

    // Constructor
    Matrix(double h_in, double dt_in, double T_in);

    // Takes indices (i, j) and returns the corresponding k index
    int index_k(int i, int j);

    // Set potential from file
    void set_potential(std::string filename);

    // Fill matrices A and B
    void fill_matrices();

    // Solve matrix equation Au^(n+1) = Bu^n
    void solve();

    // Set the initial state of the system
    void set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                           double sigma_y, double p_y);

};
