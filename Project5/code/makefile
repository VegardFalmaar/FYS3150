compile:
	g++ -o main.exe -Wall -fopenmp -O3 -std=c++17 main.cpp forward_euler.cpp backward_euler.cpp crank_nicolson.cpp solvers.cpp tridiag.cpp two_dimensions.cpp

compile_test:
	g++ -o unit_test.exe -Wall unit_test.cpp tridiag.cpp solvers.cpp forward_euler.cpp backward_euler.cpp crank_nicolson.cpp two_dimensions.cpp -fopenmp

run_test:
	./unit_test.exe

run_test_memory:
	valgrind --leak-check=yes ./unit_test.exe

test: compile_test run_test
test_memory: compile_test run_test_memory
