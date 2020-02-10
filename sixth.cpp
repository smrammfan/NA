#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

double  a = 1;
double b = 3;

double func(double x){
	return x + x*x + 4*x*log(x);
}

double* coef_Ai(int N){
	double* A = new double[N];
	double h = (b - a) / N;
	double x = a;
	for (int i = 0; i < N; ++i)	{
		A[i] = func(x);
		x += h;
	}
	return A;
}

double* coef_Di(double* C, int N){
	double* D = new double[N];
	double h = (b - a) / N;

	D[0] = C[0];
	for (int i = 1; i < N; ++i)
		D[i] = (C[i] - C[i - 1]) / h;
	return D;
}

double* coef_Bi(double* C, double* D, int N){
	double* B = new double[N];
	double h = (b - a) / N;
	B[0] = 0;

	for (int i = 1; i < N; ++i)
		B[i] = (h * C[i] / 2) - (h * h * D[i] / 2) + (func(a + h*i) - func(a + h*(i - 1))) / h;
	return B;
}

double **SLAR(int N){
	double **slar_matrix = new double*[N];
	for (int i = 0; i < N; i++)
		slar_matrix[i] = new double[N + 1];

	double h = (b - a) / N;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N + 1; j++)
			slar_matrix[i][j] = 0;
	double x_1, x_i, x_n;

	for (int i = 0; i < N; i++){
		slar_matrix[i][i] = 4 * h;
		if (i>0)
			slar_matrix[i][i - 1] = h;
		if (i < N)
			slar_matrix[i][i + 1] = h;
		x_1 = a + h * (i - 1);
		x_i = a + h * i;
		x_n = a + h * (i + 1);
		slar_matrix[i][N] = 6 * (func(x_n) - 2 * func(x_i) + func(x_1)) / h;
	}
	return slar_matrix;
}


double * Sweep_method(double **SLAR, int N){
	
	double *C = new double[N];
	double *alfa = new double[N];
	double *betta = new double[N];

	alfa[0] = -SLAR[0][1] / SLAR[0][0];
	betta[0] = SLAR[0][N] / SLAR[0][0];
	alfa[1] = -SLAR[1][2] / SLAR[1][1];
	betta[1] = SLAR[1][N] / SLAR[1][1];

	for (int i = 1; i < N; i++){
		alfa[i + 1] = -SLAR[i][i + 1] / (SLAR[i][i - 1] * alfa[i] + SLAR[i][i]);
		betta[i + 1] = (SLAR[i][N] - SLAR[i][i - 1] * betta[i]) / (SLAR[i][i - 1] * alfa[i] + SLAR[i][i]);
	}
	C[N - 1] = (SLAR[N - 1][N] - SLAR[N - 1][N - 3] * betta[N - 1]) / (SLAR[N - 1][N - 3] * alfa[N - 1] + SLAR[N - 1][N - 2]);
	for (int i = N - 1; i >= 0; i--){
		C[i] = alfa[i + 1] * C[i + 1] + betta[i + 1];
	}
	return C;
}

double Splain_Value(double x, double* A, double* B, double* C, double* D, int N){
	int i = 0;
	double h = (b - a) / N, xi = a + i * h;
	
	while (x >= xi){
		xi += h;
		i++;
	}
	i--;
	xi = a + i * h;
	return A[i] + B[i] * (x - xi) + C[i] * (x - xi) * (x - xi) / 2 + D[i] * (x - xi) * (x - xi) * (x - xi) / 6;
}

void Interpolation(int N){

	ofstream fout("lab6.txt");
	int count = (int)((b - a) * 50 + 1);
	double step = (b - a) / (double)count;
	double x, y;
	double* A, *B, *C, *D ;
	double **Matrix = new double*[N];
	for (int i = 0; i < N; i++)
		Matrix[i] = new double[N + 1];

	Matrix = SLAR(N);
	A = coef_Ai(N);
	C = Sweep_method(Matrix, N);
	D = coef_Di(C, N);
	B = coef_Bi(C, D, N);

	x = a;
	for (int i = 0; i <= count; ++i){
		y = Splain_Value(x, A, B, C, D, N);
		fout << x << ";" << y << endl;
		x += step;
	}
	fout.close();
	return;
}


int main(){
	int N = 60;
	cout << "Splain Interpolation" << endl;
	cout << "N = " << N << endl;
	Interpolation(N);
	cout << "Finished" << endl;

	return 0;
}
