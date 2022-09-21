//
// program defines a method for finding the maximum element within the off-diagonals of a symmetric matrix.
// The method is tested for a symmetric matrix A.
//
#include <armadillo>
#include <iostream>
#include <cmath>

double max_element_offdiag(const arma::mat &symmetric_matrix, int &row_max_element, int &col_max_element)
{

    const int matrix_dims = symmetric_matrix.n_cols;
    double current_element_sup;
    double max_offdiag_element = 0;
    for (int i = 0; i < matrix_dims-1; i++)
    {
        current_element_sup = symmetric_matrix(i + 1, i);
        if (std::abs(current_element_sup) > std::abs(max_offdiag_element))
        {
            row_max_element = i + 1;
            col_max_element = i;
            max_offdiag_element = symmetric_matrix(row_max_element, col_max_element);
        }
    }
    return max_offdiag_element;
};

int main()
{
    arma::mat A =
        {
            {1, 0, 0, 0.5},
            {0, 1, -0.7, 0},
            {0, -0.7, 1, 0},
            {0.5, 0, 0, 1}
        };
    std::cout << "matrix A:\n" << A <<std::endl;
    int k;
    int l;
    std::cout << "maximum off-diagonal element of A: " << max_element_offdiag(A,k,l) << ", k = " << k + 1 << ", l = " << l +1 <<  std::endl;
}