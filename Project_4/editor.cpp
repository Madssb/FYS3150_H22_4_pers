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

    arma::Mat<int> lattice;  //lagrer l x l lattice

    int l;  //l x l ising
    double M; //total magnetic
    double J;// coupling constant
    double E; //energi
    int N; //l*l
    double T; //temperatur
    double Bz; 

    int n; //antall cycles

    double avg_E, avg_EE, avg_M, avg_MM; //gjennomsnitlig energi og magnet

    void print(); //printer tilstanden av latticen med E og M

    void kjor_MCMC(); //kjører metropolis n ganger 

    void metropolis();  //utfører metorpolis på latticen


    private:

    std::map<int, double> dE;

    void Energi();

    unsigned int seed;

    mt19937 generator;

    vector<int> spin;

    void generate_lattice();

    





};

ising::ising(int n_)
{
    l = n_;

    lattice = arma::Mat<int>(l, l);
    J = 1;
    Bz = 1.380e-23; //J/K;
    M = 0.;
    E = 0.;
    N = l*l;
    n = 1000000;
    T = 1.;

    avg_E =  avg_EE =  avg_M = avg_MM =  0.;

    dE[8] = std::exp(-8 / T);
    dE[4] = std::exp(-4 / T);
    dE[0] = 1.;
    dE[-4] = 1.;
    dE[-8] = 1.;

    seed = chrono::system_clock::now().time_since_epoch().count();
  // Seed it with our seed
    generator.seed(seed);

    spin = {1, -1};
    generate_lattice();

}


void ising::metropolis()
{

    //for all spins
    std::uniform_real_distribution<double> r(0., 1.);

    for (int i = 0; i < N; i++) 
    {
        //random position
        uniform_int_distribution<int> my_01_pdf(0,l-1);
        int x = my_01_pdf(generator);
        int y = my_01_pdf(generator);

        int E1 = lattice(y, x) * lattice(y, (x+1+ l)%l) + lattice(y, x) * lattice(y, (x-1+ l)%l) + lattice(y, x) * lattice((y+1+ l)%l, x) + lattice(y, x) * lattice((y-1+ l)%l, x);

        int E2 = E1*-1;
        int deltaE = E2 - E1; 
        
        //test
        if (r(generator) <= dE[deltaE]) //accept
        {
            //update E and M 
            lattice(y, x) *= -1;
            //cout << M << "  " << (double) 2*lattice(y, x) << endl;
            M += 2*lattice(y, x);
            E += deltaE;
        }

        //else keep old konfiguaration
    }
}

void ising::generate_lattice()
{
    //unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    //mt19937 generator;
  // Seed it with our seed
    //generator.seed(seed);
    uniform_int_distribution<int> my_01_pdf(0,1);

    for (int i = 0; i < l; i++)
    {
        for (int j = 0; j < l; j++)
        {
            int val = spin[my_01_pdf(generator)];
            lattice(i, j) = val;
            M += val;
        }
    }

    Energi();

}

void ising::print()
{
    lattice.print();
    cout << "M " <<  M << endl;
    cout << "E " <<  E << endl;

}



void ising::Energi()
{ 
    for (int i = 0; i < l; i ++)
    {
        for (int j = 0; j < l; j++)
        {
            E += lattice(i, j) * lattice(i, (j+1)%l) + lattice(i, j) * lattice((i+1)%l, j);
        }
    }
}



void ising::kjor_MCMC()
{
    ofstream file("expectationvalues.txt");

    int N = l*l;
    for (int i = 1; i <= n; i++)
    {
        double norm = 1./i;
        //MCMC();
        metropolis(); //endrer E += dE og M lxl ganger 
        avg_E += E;
        avg_EE += E*E;
        avg_M += std::abs(M);
        avg_MM += M*M;

        file << i << "  " << avg_E*norm/N << "    " << avg_EE*norm/N/N << "   " << avg_M*norm/N << "    " << avg_MM*norm/N/N << endl;
        //file << i << "  " << avg_E*norm/N << "    " << avg_EE*norm/N/N << endl;

    }

    file.close();
}


int main()
{

    

    ising test_ising(2);
    test_ising.print();

    cout << "start "<< endl;

    test_ising.kjor_MCMC();

    cout << "slutt "<< endl;

    return 0;
}

//yoooooooo