all:
	g++ -o serial serial.cpp -Wall -pedantic -ansi
	g++ -o parallel parallel.cpp -fopenmp -Wall -pedantic -ansi

clean:
	rm -rf *o serial parallel
