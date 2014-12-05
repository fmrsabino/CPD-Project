all:
	mpic++ -o mpi_old mpi_old.cpp -Wall -ansi
	mpic++ -o mpi_async_old mpi_async_old.cpp -Wall -ansi
	mpic++ -o lcs-mpi lcs-mpi.cpp -Wall -ansi
	cp mpi_old /mnt/nimbus/glusterfs/cpd/3
	cp mpi_async_old /mnt/nimbus/glusterfs/cpd/3
	cp lcs-mpi /mnt/nimbus/glusterfs/cpd/3
clean:
	rm -rf *o mpi_old mpi_async lcs-mpi
