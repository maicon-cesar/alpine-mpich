//**********************************************************************************************************/
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
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <time.h>

using namespace std;

/**********************************************************************************************************/
float stof(string value)
{
    const char *s = value.c_str();
    float rez = 0, fact = 1;

    if (*s == '-') {
        s++;
        fact = -1;
    }

    for (int point_seen = 0; *s; s++){
        if (*s == '.'){
            point_seen = 1; 
            continue;
        }
    
        int d = *s - '0';
        if (d >= 0 && d <= 9) {
            if (point_seen) fact /= 10.0f;
            rez = rez * 10.0f + (float)d;
        }
    } 

    return rez * fact;
}

/**********************************************************************************************************/
int init_matrix(float *matrix, int len)
{
    string value;
    ifstream infile;

    infile.open("matrix.dat");

    if (infile.is_open()) {
        for (int i = 0; i < len; i++) {
            for (int j = 0; j < len; j++) {
                getline(infile, value);
                matrix[i * len + j] = stof(value);
            }
        }

        infile.close();
        return 0;
    }

    return 1;
}

/**********************************************************************************************************/
void print_matrix(float *matrix, int len)
{
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            cout << setprecision(3) << matrix[i * len + j] << "\t";

        }
        cout << endl;
    }
    cout << endl;
}

/**********************************************************************************************************/
void write_matrix(float *matrix, int size)
{
    ofstream outfile;
    outfile.open("matrix_mpi_gauss.dat", ios::out);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            outfile << setprecision(3) << matrix[i * size + j]  << ' ';
        }
        outfile << endl;
    }

    outfile.close();
}

/**********************************************************************************************************/
void verify_solution(float *matrix1, float *matrix2, int len)
{

    float epsilon = 0.005;

    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            assert(abs(matrix1[i * len + j] - matrix2[i * len + j]) <= epsilon);
        }
    }
}

/**********************************************************************************************************/
int main(int argc, char *argv[])
{
    double t_start;
    double t_end;
    double t_total;
    int name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int rank;
    int size;
    int res;
    int len;
    int log = 0;

    if (log)
       cout << "Iniciando..." << endl;

    if (argc < 2) {
        cout << "Uso: " << argv[0] << " len" << endl;
        return 0;
    } else {
        len = atoi(argv[1]);

        if (argc > 2)
            log = 1;
    }

    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(processor_name, &name_len);

    if (log) {
        cout << endl << "Hello do processador " << processor_name << " rank " << rank << 
            " de " << size  << " processadores." << endl;
    }

    int num_rows = len / size;  // número de linhas baseado no número de ranks.

    float *matrix;

    if (rank == 0) {
        srand(time(NULL));
 
        matrix = new float [len * len];

        res = init_matrix(matrix, len);
        if (res) {
            MPI_Finalize();
            delete [] matrix;
            exit(0);
        }
    }
    
    float *sub_matrix = new float[len * num_rows];

    MPI_Scatter(matrix, len * num_rows, MPI_FLOAT, sub_matrix,
            len * num_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float *row = new float[len];

    if (rank == 0) {
        t_start = MPI_Wtime();
    }

    int pivot;
    int scale;
    int column;
    int start_row;

    // Receptores aqui
    start_row = rank * num_rows;
    for (int i = 0; i < start_row; i++) {
        // Aguarda uma linha
        MPI_Bcast(row, len, MPI_FLOAT, i / num_rows, MPI_COMM_WORLD);

        for (int j = 0; j < num_rows; j++) {
            scale = sub_matrix[j * len + i];

            for (int k = i + 1; k < len; k++) {
                sub_matrix[j * len + k] -= scale * row[k];
            }

            sub_matrix[j * len + i] = 0;
        }
    }

    // Envio daqui
    for (int i = 0; i < num_rows; i++) {

        column = rank * num_rows + i;
        pivot = sub_matrix[i * len + column];

        for (int j = column + 1; j < len; j++) {
           sub_matrix[i * len + j] /= pivot;
        }

        sub_matrix[i * len + column] = 1;

        // Preenche linha a ser enviada
        memcpy(row, &sub_matrix[i * len], len * sizeof(float));

        // Broadcast linha normalizada
        MPI_Bcast(row, len, MPI_FLOAT, rank, MPI_COMM_WORLD);

        // Atualiza o resto das linhas nesse rank
        for (int j = i + 1; j < num_rows; j++) {
            scale = sub_matrix[j * len + column];

            for (int k = column + 1; k < len; k++) {
                sub_matrix[j * len + k] -= scale * row[k];
            }
            
            sub_matrix[j * len + column] = 0;
        }
    }

    // Sincroniza finalização do rank
    for (int i = rank * num_rows + 1; i < len; i++) {
        MPI_Bcast(row, len, MPI_FLOAT, i / num_rows, MPI_COMM_WORLD);
    }

    // Barrier para esperar fim dos cálculos
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        t_end = MPI_Wtime();
        t_total = t_end - t_start;
    }

    // Coleta todas as sub-matrizes
    MPI_Gather(sub_matrix, len * num_rows, MPI_FLOAT, matrix,
            len * num_rows, MPI_FLOAT, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    if (rank == 0) {
        //print_matrix(matrix, len);
        cout << t_total << " segundos." << endl;
        write_matrix(matrix, len);
    }

    if (rank == 0) {
        delete[] matrix;
    }
    delete[] sub_matrix;
    delete[] row;

    return 0;
}
