all: mpi_cosopt
mpi_cosopt: mpi_cosopt.c
	/home3/ishiwata/mpich2/bin/mpicc mpi_cosopt.c -o mpi_cosopt -lm -ldsfmt -lsfmt -L/home3/ishiwata/lib -L/home3/ishiwata/mpich2/lib -I/home3/ishiwata/include -I/home3/ishiwata/mpich2/include -O3
clean: mpi_cosopt
	rm mpi_cosopt
