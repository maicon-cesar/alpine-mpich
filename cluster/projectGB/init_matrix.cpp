#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>


using namespace std;

/**********************************************************************************************************/
void init_matrix(float *matrix, int size)
{
    ofstream outfile;
    outfile.open("matrix.dat", ios::out);

    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            matrix[i * size + j] = (float(rand())/float(RAND_MAX)) * (100 - -100) + -100;

            outfile << setprecision(3) << matrix[i * size + j]  << endl;
        }
    }

    outfile.close();
}

/**********************************************************************************************************/
int main(int argc, char *argv[])
{
    float *matrix;
    int size;

    if (argc != 2) {
        cout << "Usage: " << argv[0] << " len" << endl;
        return 0;
    }

    try {
        size = atoi(argv[1]);
    } catch (...) {
        cout << "Tamanho de matriz informado invalido!" << endl;
        return 0;
    }

    matrix = new float [size * size];
    init_matrix(matrix, size);

    delete[] matrix;
    return 0;
}
