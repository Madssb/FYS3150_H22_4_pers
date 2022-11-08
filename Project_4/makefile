all: compile run

compile:
	g++ -O2 ising.cpp src/lattice.cpp src/mcmc.cpp src/energy_magnetization.cpp src/expectation_values.cpp -std=c++11 -I include -o ising.exe -larmadillo

run:
	./ising.exe
