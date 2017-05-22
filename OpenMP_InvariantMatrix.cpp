#include <QCoreApplication>
#include <QDebug>
#include <QTime>
#include <iostream>
#include <omp.h>

#define THREAD_NUMS 4

void inversion(double **A, int N)
{
    double temp;
    double **E = new double *[N]; // E - matrix

    for (int i = 0; i < N; i++)
        E[i] = new double [N];
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            E[i][j] = 0.0;
            if (i == j)
                E[i][j] = 1.0;
        }
    for (int k = 0; k < N; k++)
    {
        temp = A[k][k];
        #pragma omp parallel for
        for (int j = 0; j < N; j++)
        {
            A[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < N; i++)
        {
            temp = A[i][k];
            #pragma omp parallel for
            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }
    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = A[i][k];
            #pragma omp parallel for
            for (int j = 0; j < N; j++)
            {
                A[i][j] -= A[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = E[i][j];

    for (int i = 0; i < N; i++)
        delete [] E[i];

    delete [] E;
}

void GetMatr(double **mas, double **p, int i, int j, int m) // ex i string j column
{
  int ki, kj, di, dj;
  di = 0;
  for (ki = 0; ki<m - 1; ki++) { // проверка индекса строки
    if (ki == i) di = 1;
    dj = 0;
    for (kj = 0; kj<m - 1; kj++) { // проверка индекса столбца
      if (kj == j) dj = 1;
      p[ki][kj] = mas[ki + di][kj + dj];
    }
  }
}

void printMatrix(double **matrix, int size) {
  for (int i = 0; i<size; i++) {
    for (int j = 0; j<size; j++)
      std::cout << matrix[i][j] << " ";
    std::cout << std::endl;
  }
}

int determinant(double **mas, int m) {
  int i, det, k, n;
  double **p;
  p = new double*[m];
  for (i = 0; i<m; i++)
    p[i] = new double[m];
  det = 0;
  k = 1; //(-1) в степени i
  n = m - 1;
  if (m<1)
  {
      return 0;
  }
  if (m == 1) {
    det = mas[0][0];
    return det;
  }
  if (m == 2) {
    det = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]);
    return det;
  }
  if (m>2) {
    #pragma omp parallel for
    for (i = 0; i<m; i++) {
      GetMatr(mas, p, i, 0, m);
      det = det + k * mas[i][0] * determinant(p, n);
      k = -k;
    }
  }
  return det;
}

int main()
{
    double start_time, end_time, tick;
    int N;
    std::cout << "Enter N: ";
    std::cin >> N;
    //srand(121);
    omp_set_num_threads(THREAD_NUMS);
    start_time = omp_get_wtime();
    double **matrix = new double *[N];

    for (int i = 0; i < N; i++)
        matrix[i] = new double [N];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            matrix[i][j] = rand() % 10;
        }
    printMatrix(matrix, N);
    if(determinant(matrix, N) != 0)
        inversion(matrix, N);
    printMatrix(matrix, N);
    for (int i = 0; i < N; i++)
        delete [] matrix[i];
    delete [] matrix;
       end_time = omp_get_wtime();
    qDebug("Time elapsed: %lf\n", end_time-start_time);
    std::cin.get();
    return 0;
}
