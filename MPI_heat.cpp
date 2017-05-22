#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

#define NUMBER_OF_NODES  300000

void Initialize(double* T0, double* T1, int count)
{
	for (int i = 0; i < count; ++i)
	{
		T0[i] = 100;
		T1[i] = 100;
	}
}
void Swap(double** T0, double** T1)
{
	double* tmp;
	tmp = (*T0);
	(*T0) = (*T1);
	(*T1) = tmp;
}

void Step(const double* T0, double* T1, const int& size, const double& rTLeft, const double& rTRight, const double& dx,
	const double& dt, const int& r, const int& pc)
{
	//std::cout << "Step My rTLeft is " << rTLeft << " " << rTRight << std::endl;
	for (int i = 0; i < size; ++i)
	{
		if (i == 0)
		{
			if (r == 0)
				T1[i] = rTLeft;
			else
				T1[i] = T0[i] + dt / (dx*dx) * (rTLeft - 2 * T0[i] + T0[i + 1]);
			//std::cout << T1[i] << std::endl;
			//T1[i] = rTLeft;
		}
		else if (i == size - 1)
			if (r == pc - 1)
				T1[i] = rTRight;
			else
				T1[i] = T0[i] + dt / (dx*dx) * (T0[i - 1] - 2 * T0[i] + rTRight);
			
//			T1[i] = rTRight;
		else
			T1[i] = T0[i] + dt / (dx*dx) * (T0[i - 1] - 2 * T0[i] + T0[i + 1]);
	}
	/*if (r == 1)
		for (int i = 0; i < size - 1; ++i)
			std::cout << T1[i] << " ";
	std::cout << std::endl;*/
}
void Send(const double& sTLeft, const double& sTRight, double& rTLeft, double& rTRight, const int& r, const int& pc)
{
	//if(r == 0)
	//std::cout << "Send My rank is " << r;
	MPI_Status status;

	if (r < pc - 1)
		MPI_Send(reinterpret_cast<void*>(const_cast<double*>(&sTRight)), 1, MPI_DOUBLE, r + 1, 1, MPI_COMM_WORLD);
	if (r > 0)
		MPI_Recv(reinterpret_cast<void*>(&rTLeft), 1, MPI_DOUBLE, r - 1, 1, MPI_COMM_WORLD, &status);

	if (r > 0)
		MPI_Send(reinterpret_cast<void*>(const_cast<double*>(&sTLeft)), 1, MPI_DOUBLE, r - 1, 1, MPI_COMM_WORLD);
	if (r < pc - 1)
		MPI_Recv(reinterpret_cast<void*>(&rTRight), 1, MPI_DOUBLE, r + 1, 1, MPI_COMM_WORLD, &status);
}
void Output(double* T0, const int& size, const int& rank)
{
	char name[10];
	sprintf(name, "rank%i", rank);
	//std::cout << "My name is " << name << std::endl;
	FILE* file = fopen(name, "a");
	for (int i = 0; i < size; ++i)
	{
		if(i != size - 1)
			fprintf(file, "%lf\t", T0[i]);
		else
			fprintf(file, "%lf\t\n", T0[i]);
	}
	fclose(file);
}
void Calculate(const double* TInit, const double& a, const double& b, const double& N, const int& r, const int& pc)
{
	// Nr - node count
	//std::cout << "Calculate My rank = " << r << std::endl;
	int Nr = N / pc;
	double dx = (b - a) / N;//(Nr - 1);
	std::cout << "dx = " << dx << std::endl;
	double dt = 0.00001;
	int i1 = Nr * r;
	std::cout << "My i1 = " << i1 << std::endl;
	if (r == pc - 1)
		Nr = N - i1;
	int i2 = i1 + Nr;
	std::cout << "My i2 = " << i2 << std::endl;
	int size = Nr;
	double *T0 = new double[size];
	double *T1 = new double[size];
	//std::cout << TInit[6] << std::endl;
	for (int i = i1; i < i2; ++i)
	{
		T0[i - i1] = TInit[i1*0 + i];
		T1[i - i1] = TInit[i1*0 + i];
		//std::cout << "My T0 = " << T0[i - i1] << " Tinit = " << TInit[i1*0 + i]  << " " << i << std::endl;
	}	
	for (int s = 0; s < 2500; ++s)
	{			
		//std::cout << "Cicle " << r << std::endl;
		double sTLeft = T0[0], sTRight = T0[size - 1];
		double rTLeft = T0[0], rTRight = T0[size - 1];
		Send(sTLeft, sTRight, rTLeft, rTRight, r, pc);
		//std::cout << "Rank = " << r << "STL "<<sTLeft << " " << sTRight << " " << sTLeft << " " << rTRight << std::endl;
		//std::cout << "Send OK" << std::endl;
		Step(T0, T1, size, rTLeft, rTRight, dx, dt, r, pc);
		if(s % 10 == 0)
			//Output(T1, size, r);
		Swap(&T0, &T1);
		;
	}
	std::cout << "Finish " << r <<std::endl;
}

int main(int argc, char** argv)
{
	float fTimeStart = clock() / (float)CLOCKS_PER_SEC;
	double a = 0;
	double b = 1;

	MPI_Init(&argc, &argv);
	int process_count = 0;
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &process_count);

	double *TInit = new double[NUMBER_OF_NODES];
	for (int i = 1; i < NUMBER_OF_NODES - 1; ++i)
	TInit[i] = 0;
	TInit[0] = 0;
	TInit[NUMBER_OF_NODES - 1] = 100;
//	Output(TInit, NUMBER_OF_NODES, -100);
//	std::cout << "Hi" << std::endl;
	Calculate(TInit, a, b, NUMBER_OF_NODES, rank, process_count);

	
	MPI_Finalize();
	float fTimeStop = clock() / (float)CLOCKS_PER_SEC;
	std::cout << "Working time = " << fTimeStop-fTimeStart<< std::endl;
	delete[] TInit;
	return 0;
}
