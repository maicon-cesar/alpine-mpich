/**********************************************************************************************************/
// Implementa eliminação Gaussiana usando programação paralela (MPI)
// Projeto GB - Arquitetura de Computadores II
// Professor Lúcio Rene Prade
// Maicon Cesar Canales de Oliveira
// Filipe Schenkel de Souza
// Vitor de Castro
/**********************************************************************************************************/

#include <stdlib.h>
#include <mpi.h>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <assert.h>
#include <time.h>

using namespace std;

/**********************************************************************************************************/
// Função que implementa de modo serial para comparação
/**********************************************************************************************************/
void get_serial(float *matrix, int n){

    float pivot;

    for (int i = 0; i < n - 1; i++) {
        pivot = matrix[i * n + i];

        for (int j = i + 1; j < n; j++) {
            matrix[i * n + j] /= pivot;
        }

        matrix[i * n + i] = 1.0;

        float scale;
        for (int j = i + 1; j < n; j++) {
            scale = matrix[j * n + i];

            for (int k = i + 1; k < n; k++) {
                matrix[j * n + k] -= matrix[i * n + k] * scale;
            }

            matrix[j * n + i] = 0;
        }
    }

    matrix[(n - 1) * n + n - 1] = 1;
}

/**********************************************************************************************************/
void init_matrix(float *matrix, int N)
{
    srand(time(NULL));
    
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            matrix[i * N + j] = (float(rand())/float(RAND_MAX)) * (100 - -100) + -100; 
        }
    }
}

/**********************************************************************************************************/
void print_matrix(float *matrix, int N)
{
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            cout << setprecision(3) << matrix[i * N + j] << "\t";

        }
        cout << endl;
    }
    cout << endl;
}

/**********************************************************************************************************/
void verify_solution(float *matrix1, float *matrix2, int N)
{

    float epsilon = 0.005;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            assert(abs(matrix1[i * N + j] - matrix2[i * N + j]) <= epsilon);
        }
    }
}

/**********************************************************************************************************/
int main(int argc, char *argv[])
{
    int N = 1024;

    double t_start;
    double t_end;
    double t_total;
    int name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int rank;
    int size;

    cout << "Iniciando..." << endl;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(processor_name, &name_len);

    cout << endl << "Hello do processador " << processor_name << " rank " << rank << 
    " de " << size  << " processadores." << endl;

    int num_rows = N / size;  // número de linhas baseado no número de ranks.

    float *matrix;

    if (rank == 0) {
        matrix = new float [N * N];

        init_matrix(matrix, N);
    }
    
    float *sub_matrix = new float[N * num_rows];

    MPI_Scatter(matrix, N * num_rows, MPI_FLOAT, sub_matrix,
            N * num_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float *row = new float[N];

    if (rank == 0) {
        t_start = MPI_Wtime();
    }

    int pivot;
    int scale;
    int column;
    int start_row;

    // Receptores aqui
    start_row = rank * num_rows;
    for(int i = 0; i < start_row; i++){
        // Wait for the preceeding ranks to forward us a row

        // Aguarda uma linha
        MPI_Bcast(row, N, MPI_FLOAT, i / num_rows, MPI_COMM_WORLD);

        for(int j = 0; j < num_rows; j++) {
            scale = sub_matrix[j * N + i];

            for(int k = i + 1; k < N; k++) {
                sub_matrix[j * N + k] -= scale * row[k];
            }

            sub_matrix[j * N + i] = 0;
        }
    }

    // Envio daqui
    for (int i = 0; i < num_rows; i++) {

        column = rank * num_rows + i;
        pivot = sub_matrix[i * N + column];

        for (int j = column + 1; j < N; j++) {
           sub_matrix[i * N + j] /= pivot;
        }

        sub_matrix[i * N + column] = 1;

        // Preenche linha a ser enviada
        memcpy(row, &sub_matrix[i * N], N * sizeof(float));

        // Broadcast linha normalizada
        MPI_Bcast(row, N, MPI_FLOAT, rank, MPI_COMM_WORLD);

        // Atualiza o resto das linhas nesse rank
        for(int j = i + 1; j < num_rows; j++){
            scale = sub_matrix[j * N + column];

            for(int k = column + 1; k < N; k++){
                sub_matrix[j * N + k] -= scale * row[k];
            }
            
            sub_matrix[j * N + column] = 0;
        }
    }

    // Sincroniza finalização do rank
    for (int i = rank * num_rows + 1; i < N; i++) {
        MPI_Bcast(row, N, MPI_FLOAT, i / num_rows, MPI_COMM_WORLD);
    }

    // Barrier para esperar fim dos cálculos
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        t_end = MPI_Wtime();
        t_total = t_end - t_start;
    }


    // Coleta todas as sub-matrizes
    MPI_Gather(sub_matrix, N * num_rows, MPI_FLOAT, matrix,
            N * num_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    if (rank == 0) {
        //print_matrix(matrix, N);
        cout << t_total << " segundos." << endl;
    }

    if (rank == 0) {
        delete[] matrix;
    }
    delete[] sub_matrix;
    delete[] row;

    return 0;
}
