#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
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
void init_matrix(float *matrix, int len)
{
    srand(time(NULL));

    for(int i = 0; i < len; i++) {
        for(int j = 0; j < len; j++) {
            matrix[i * len + j] = (float(rand())/float(RAND_MAX)) * (100 - -100) + -100;
        }
    }
}

/**********************************************************************************************************/
void read_matrix(float *matrix, int len)
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
    }
}

/**********************************************************************************************************/
void write_matrix(float *matrix, int size)
{
    ofstream outfile;
    outfile.open("matrix_gauss.dat", ios::out);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            outfile << setprecision(3) << matrix[i * size + j]  << ' ';
        }
        outfile << endl;
    }

    outfile.close();
}

/**********************************************************************************************************/
void gauss_serial(float *matrix, int len){

    float pivot;

    for (int i = 0; i < len - 1; i++) {
        pivot = matrix[i * len + i];

        for (int j = i + 1; j < len; j++) {
            matrix[i * len + j] /= pivot;
        }

        matrix[i * len + i] = 1.0;

        float scale;
        for (int j = i + 1; j < len; j++) {
            scale = matrix[j * len + i];

            for (int k = i + 1; k < len; k++) {
                matrix[j * len + k] -= matrix[i * len + k] * scale;
            }

            matrix[j * len + i] = 0;
        }
    }

    matrix[(len - 1) * len + len - 1] = 1;
}

/**********************************************************************************************************/
void print_matrix(float *matrix, int len)
{
    for(int i = 0; i < len; i++) {
        for(int j = 0; j < len; j++) {
            cout << setprecision(3) << matrix[i * len + j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

/**********************************************************************************************************/
int main(int argc, char *argv[])
{
    int len;
    float *matrix;

    if (argc != 2) {
        cout << "Uso: " << argv[0] << " len" << endl;
        return 0;
    } else {
        len = atoi(argv[1]);
    }

    const clock_t begin_time = clock();

    matrix = new float[len * len];
    read_matrix(matrix, len);

    gauss_serial(matrix, len);

    write_matrix(matrix, len);

    delete[] matrix;

    cout << float(clock() - begin_time) / CLOCKS_PER_SEC << " segundos." << endl;    
    return 0;
}
