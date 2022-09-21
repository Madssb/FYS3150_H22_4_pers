#include <armadillo>
#include <iostream>

double max_offdiag_symmetric(const arma::mat &A, int &k, int &l)
{
    int matrix_dims = A.n_rows();
    double current_element_sup;
    double max_element_sup = 0; 
    for (indice = 0; i < matrix_dims() - ; i++)
    {
        current_element_sup = A(i + 1, i)
        if(current_element_sup > max_element_sup)
        {
            row_largest = i+1;
            col_largest = i;
            max_element_sup = A(i+1,i);
        }
    }
    k = row_largest;
    l = col_largest;
}