/*
This file contains the definition of the class Matrix
*/

class Matrix
{
  public:

    // Constructor
    Matrix(const int M_in);

    // Creating a matrix
    arma::mat create_matrix();

    // Takes indices (i, j) and returns the corresponding k index
    int index_k(int i, int j);

    //
};
