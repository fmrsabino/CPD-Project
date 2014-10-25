all:
	g++ -o serial serial.cpp -Wall -pedantic -ansi
	g++ -o parallel parallel.cpp -Wall -pedantic -ansi
run:
	./serial ex48k.30k.in
clean:
	rm -rf *o serial parallel
