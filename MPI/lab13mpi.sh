for i in 1 2 3 4 5 6 7 8 9 10
do
	scp mpi_fix ist173239@lab13p$i.rnl.ist.utl.pt:/tmp
	scp my-nodes.txt ist173239@lab13p$i.rnl.ist.utl.pt:/tmp
	scp -r ./test-files/ ist173239@lab13p$i.rnl.ist.utl.pt:/tmp
done