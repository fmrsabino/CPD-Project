all:
	g++ -o serial serial.cpp -Wall -pedantic -ansi
	g++ -o parallel parallel.cpp -fopenmp -Wall -pedantic -ansi
runs:
	./serial ex3k.8k.in

runp:
	./parallel ex3k.8k.in
clean:
	rm -rf *o serial parallel
