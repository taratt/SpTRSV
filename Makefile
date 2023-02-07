CC=g++-12
CFLAGS=-std=c++11 -fopenmp


# It is assumed that the boost library is installed and is in the /usr/local/include/boost/ path
main: main.cpp inputFormatter.h DependencyGraphCalulations.h matrixCalculations.h outputVerifier.h
	$(CC) $(CFLAGS) main.cpp -o main -I /usr/local/include/boost/
