#include <iostream>
#include <fstream>
using namespace std;
#include <cmath>
#include <cstdlib> 

#define ARMA_DONT_USE_STD_MUTEX

#include <armadillo>



int main()
{
    int A = 20;

    int* ptr; 
    ptr = &A;
    cout << ptr << " " << *ptr << endl;
    cout << A  << endl;

    cout << "change" << endl;
    *ptr = 10;

    cout << ptr << " " << *ptr << endl;
    cout << A  << endl;
    
    return 0;
}