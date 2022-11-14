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

    arma::Mat<int> lattice;  //lagrer N x N lattice

    int l;  //l x l ising
    int M; //total magnetic
    double J;// coupling constant
    double E; //energi
    int N; //l*l
    double T; //temperatur

    int n; //antall cycles

    double avg_E, avg_EE, avg_M, avg_MM; 



    double Energi_; //total energi fra alle flippene 

    double Bz; 

    void print();
    void flip();  //flipper en tilfeldig spin

    

    void MCMC();
    void kjor_MCMC();

    int MCMC_index;

    int flip_index;
    vector<double> acceptet_states;
    vector<double> Es;  // energiene fra markov-kjeden
    vector<int> Bs; //total magnet fra markov-kjeden


    void expectationvalue();

    void metropolis();
    private:

    std::map<int, double> dE;


    void Energi();
    void Energi2x2();

    int Energi2x2_();

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
    n = 700;
    T = 1.;

    avg_E =  avg_EE =  avg_M = avg_MM =  0;

    dE[8] = std::exp(-8 / T);
    dE[4] = std::exp(-4 / T);
    dE[0] = 1.;
    dE[-4] = 1.;
    dE[-8] = 1.;


    MCMC_index = 0;
    

    seed = chrono::system_clock::now().time_since_epoch().count();
  // Seed it with our seed
    generator.seed(seed);

    spin = {1, -1};
    generate_lattice();

    Es = vector<double>(n, 0.);
    Bs = vector<int>(n, 0);
    //acceptet_states = vector<int>(n*N*N, 0);
    flip_index = 0;
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
        //cout << "y, x: " <<  y << " " << x << endl;

        //energy difference 
        //cout << lattice(y, x) << " " << lattice(y, (x+1 + l)%l)  << endl;
        //cout << lattice(y, x) << " " << lattice(y, (x-1 + l)%l) << endl;
        //cout << lattice(y, x) << " " << lattice((y+1+ l)%l, x) << endl;
        //cout << lattice(y, x) << " " << lattice((y-1+ l)%l, x) << endl;

        //int E1 = Energi2x2_();
        int E1 = lattice(y, x) * lattice(y, (x+1+ l)%l) + lattice(y, x) * lattice(y, (x-1+ l)%l) + lattice(y, x) * lattice((y+1+ l)%l, x) + lattice(y, x) * lattice((y-1+ l)%l, x);

        int E2 = E1*-1;
        int deltaE = E2 - E1; 
        //cout << deltaE << endl;
        
        //test
        if (r(generator) <= dE[deltaE]) //accept
        {
            //update E and M 
            lattice(y, x) *= -1;
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

    //Bs[MCMC_index] = M;
    Energi();
    //Energi2x2();
}

void ising::print()
{
    lattice.print();
    cout << "M " <<  M << endl;
    cout << "E " <<  E << endl;

}

void ising::flip()
{

    uniform_int_distribution<int> my_01_pdf(0,l-1);
    int nr1 = my_01_pdf(generator);
    int nr2 = my_01_pdf(generator);

    lattice(nr1,nr2 ) *= -1;
    M += 2*lattice(nr1,nr2 );

    Energi2x2();

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

void ising::Energi2x2()
{
    E = lattice(0, 0)* lattice(0, 1) + lattice(0, 0)* lattice(1, 0) + lattice(1, 1)* lattice(0, 1) + lattice(1, 1)* lattice(1, 0);
}

int ising::Energi2x2_()
{
    return lattice(0, 0)* lattice(0, 1) + lattice(0, 0)* lattice(1, 0) + lattice(1, 1)* lattice(0, 1) + lattice(1, 1)* lattice(1, 0);
}


void ising::kjor_MCMC()
{
    //ofstream file("expectationvalues.txt");

    int N = l*l;
    for (int i = 1; i < n; i++)
    {
        double norm = 1./i;
        //MCMC();
        metropolis(); //endrer E += dE lxl ganger 
        avg_E += E;
        avg_EE += E*E;
        avg_M += M;
        avg_M += M*M;

        cout << i << "  " << avg_E*norm/N << "    " << avg_EE*norm/N/N << "   " << avg_M*norm/N/N << "    " << avg_MM*norm/N/N << endl;

    }

    //file.close();
}

void ising::expectationvalue()
{
    ofstream file("expectationvalues.txt");

    double E_ = 0. ;
    for (int i = 1; i < MCMC_index ; i++)
    {
        E_ += Es[i];

        file << i << " " << E_/i/(l*l) << endl;
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