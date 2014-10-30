all:
	g++ -o lcs-serial lcs-serial.cpp -Wall -pedantic -ansi
	g++ -o lcs-omp lcs-omp.cpp -fopenmp -Wall -pedantic -ansi
	g++ -o lcs-omp-functions lcs-omp-functions.cpp -fopenmp -Wall -pedantic -ansi

clean:
	rm -rf *o lcs-serial lcs-omp lcs-omp-functions
