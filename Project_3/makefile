all_single: compile_single run_single

compile_single:
	g++ testing_single.cpp src/particle.cpp src/penning_trap.cpp src/analytical.cpp -std=c++11 -I include -o testing_single.exe -larmadillo

run_single:
	./testing_single.exe > data_test_single.txt

all_double: compile_double run_double

compile_double:
	g++ testing_double.cpp src/particle.cpp src/penning_trap.cpp src/analytical.cpp -std=c++11 -I include -o testing_double.exe -larmadillo

run_double:
	./testing_double.exe > data_test_double.txt

all_relerr: compile_relerr run_relerr

compile_relerr:
	g++ rel_error.cpp src/particle.cpp src/penning_trap.cpp src/analytical.cpp -std=c++11 -I include -o rel_error.exe -larmadillo

run_relerr:
	./rel_error.exe

all: compile run

compile:
	g++ main.cpp src/particle.cpp src/penning_trap_2.cpp -std=c++11 -I include -o main.exe -larmadillo

run:
	./main.exe
