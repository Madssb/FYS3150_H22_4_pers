#include <armadillo>
#include <iostream>
#include <cmath>

double maximum_element(const arma::mat &symmetric_matrix, int &k, int &l)
{
    int matrix_dims = symmetric_matrix.n_rows();
    double current_element_sup;
    double max_element_sup = 0; 
    for (i = 0; i < symmetric_matrix.n_cols()  ; i++)
    {
        current_element_sup = symmetric_matrix(i + 1, i)
        if(std::abs(current_element_sup) > std::abs(max_element_sup))
        {
            //k and l are decided by the position of the largest element
            k = i+1;
            l = i;
            max_element_sup = symmetric_matrix(i+1,i);
        }
    }
    return max_element_sup
}
