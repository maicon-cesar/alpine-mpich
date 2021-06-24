all: app_mpi app init_matrix

app_mpi:
	mpic++ -o app_mpi_gauss app_mpi_gauss.cpp

app:
	mpic++ -o app_gauss app_gauss.cpp

init_matrix:
	mpic++ -o init_matrix init_matrix.cpp

clean:
	rm -rf *.o *~ app_mpi_gauss app_gauss init_matrix
