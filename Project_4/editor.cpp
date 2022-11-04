#include <iostream>
#include <fstream>
using namespace std;
#include <cmath>
#include <cstdlib> 

#define ARMA_DONT_USE_STD_MUTEX

#include <armadillo>
using namespace arma;

class ising
{
    public:


    ising(int n);

    arma::Mat<int> lattice;

    int N;  //N x N ising
    int M; //total magnetic
    double J;// coupling constant
    double E; //energi


    void print();
    void flip();

    private:


    unsigned int seed;
    mt19937 generator;



    vector<int> spin;
    void generate_lattice();





};

ising::ising(int n)
{
    N = n;

    lattice = arma::Mat<int>(N, N);
    J = 1;
    M = 0;
    E = 0;

    seed = chrono::system_clock::now().time_since_epoch().count();
  // Seed it with our seed
    generator.seed(seed);


    spin = {1, -1};
    generate_lattice();
}

void ising::generate_lattice()
{
    //unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    //mt19937 generator;
  // Seed it with our seed
    //generator.seed(seed);
    uniform_int_distribution<int> my_01_pdf(0,1);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            int val = spin[my_01_pdf(generator)];
            lattice(i, j) = val;
            M += val;

        }
    }

}

void ising::print()
{
    lattice.print();
    cout << "M " <<  M << endl;
    cout << "E " <<  E << endl;

}

void ising::flip()
{

    uniform_int_distribution<int> my_01_pdf(0,N-1);
    int nr1 = my_01_pdf(generator);
    int nr2 = my_01_pdf(generator);

    lattice(nr1,nr2 ) *= -1;
    M += 2*lattice(nr1,nr2 );


}

int main()
{

    ising test_ising(2);

    test_ising.print();
    test_ising.flip();
    test_ising.print();

    return 0;
}

