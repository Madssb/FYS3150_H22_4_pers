build: compile link

compile:
	g++ -c -O2 ising.cpp src/lattice.cpp src/mcmc.cpp -I include -std=c++11

link:
	g++ ising.o lattice.o mcmc.o -o ising.exe -larmadillo

build-parallel:
	g++ -O2 ising.cpp src/lattice.cpp src/mcmc.cpp -I include -std=c++11 -Xpreprocessor -fopenmp -o ising.exe -larmadillo -lomp
