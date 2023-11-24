#!/bin/sh

mpirun -n $((nproc/2)) bil $@ -mpi_linear_solver_server
#mpirun -n 4 bil $@ -mpi_linear_solver_server

