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

    int N;  //N x N ising
    int M; //total magnetic
    double J;// coupling constant
    double E; //energi

    int n; //antall skjeder

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


    private:


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
    N = n_;

    lattice = arma::Mat<int>(N, N);
    J = 1;
    Bz = 1.380e-23; //J/K;
    M = 0.;
    E = 0.;
    n = 10000;

    MCMC_index = 0;
    

    seed = chrono::system_clock::now().time_since_epoch().count();
  // Seed it with our seed
    generator.seed(seed);

    spin = {1, -1};
    generate_lattice();

    Es = vector<double>(n, 0.);
    Bs = vector<int>(n, 0);
    //acceptet_states = vector<int>(n*N, 0);
    flip_index = 0;
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

    //Bs[MCMC_index] = M;


    //Energi();
    Energi2x2();
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

    Energi2x2();

}

void ising::Energi()
{

    for (int i = 0; i < N; i ++)
    {
        for (int j = 0; j < N; j++)
        {
            E = lattice(i, j) * lattice(i, (j+1)%N) + lattice(i, j) * lattice((i+1)%N, j);
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


void ising::MCMC()
{

    int energi = 0;
    for (int a = 0; a < N ; a++)
    {

    //the intial state is already random

    //flip it 
    uniform_int_distribution<int> my_01_pdf(0,N-1);
    int i = my_01_pdf(generator);
    int j = my_01_pdf(generator);

    //for N
    //int E1 = lattice(i, j) * lattice(i, (j+1)%N) + lattice(i, j) * lattice(i, (j-1)%N) + lattice(i, j) * lattice((i+1)%N, j) + lattice(i, j) * lattice((i-1)%N, j);
    //for N = 2
    int E1 = Energi2x2_();

    lattice(i,j ) *= -1;
    //M += 2*lattice(i,j );

    int E2 = Energi2x2_();

    int dE = (E2 - E1);
    //E += dE;

    //latticen er flippa og vi har styr over energi forskjellen 

    double p_diff = exp(-Bz*dE);

    uniform_real_distribution<double> r(0,1);


    double A;
    if (p_diff > 1) {A = 1;}
    else {A = p_diff;}

    if (r(generator) <= A  )
    {
        //this is the new state
        //lattice(i,j ) *= -1;
        M += 2*lattice(i,j );

        //int dE = (E2 - E1);
        E = E2; // E + dE;
        energi += E2;
        acceptet_states.push_back(E*1.);
        flip_index += 1;
        //E += dE;

        //Es[MCMC_index] = E;
        //Bs[MCMC_index] = M;
        //MCMC_index += 1;
    }

    else
    {
        lattice(i,j ) *= -1;
        E = Energi2x2_();
        energi += Energi2x2_();

        acceptet_states.push_back(E*1.);
        flip_index += 1;
    }

    }
    Es[MCMC_index] = energi; ;//energi/N;
    Bs[MCMC_index] = M;
    MCMC_index += 1;

    //cout << E << endl;




}

void ising::kjor_MCMC()
{
    for (int i = 0; i < n; i++)
    {
        MCMC();
    }
}

void ising::expectationvalue()
{
    ofstream file("expectationvalues.txt");

    

    for (int i = 1; i < MCMC_index; i++)
    {
        double energi = 0;
        for (int j = 0; j < i; j++)
        {
            energi += Es[j];

        }

        file << i << " " << energi/i << endl;
    }
    
    
    file.close();

}

int main()
{

    cout << "start "<< endl;

    ising test_ising(2);
    test_ising.print();
    cout << "middel "<< endl;
    test_ising.kjor_MCMC();


    test_ising.expectationvalue();
    
    cout << test_ising.acceptet_states[2]; 

    ofstream file("Es.txt");
    for (double data : test_ising.acceptet_states)
    {
        file << data << endl;
    }
    file.close();




    cout << "slutt "<< endl;

    return 0;
}

//yoooooooo