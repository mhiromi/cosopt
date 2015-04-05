#!/bin/sh

export LD_LIBRARY_PATH=/home3/ishiwata/mpich2/lib:/home3/ishiwata/lib:$LD_LIBRARY_PATH
#$ -pe mpich 10
/home3/ishiwata/mpich2/bin/mpirun -n $NSLOTS ./mpi_cosopt $1 $2
