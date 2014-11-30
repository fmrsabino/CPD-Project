all:
	mpic++ -o mpi mpi.cpp -Wall -ansi
	mpic++ -o mpi_async mpi_async.cpp -Wall -ansi
	mpic++ -o mpi_fix mpi_fix.cpp -Wall -ansi
	cp mpi /mnt/nimbus/glusterfs/cpd/3
	cp mpi_async /mnt/nimbus/glusterfs/cpd/3
clean:
	rm -rf *o mpi
