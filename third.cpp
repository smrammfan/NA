#include <iostream>
#include <cmath>

#define N 4

double arr[N][N + 1] = { 15, 8, 19, 7, 69,
						15, 41, 16, 9, 76,
						0, 3, 20, 16, 100,
						15, 7, 8, 18, 113};

double modern_arr[N][N + 1] = {30, 12, 7, 9, 82,
							   15, 41, 16, 9, 76,
							   0, 3, 20, 16, 100,
							   0, -2, -9, -27, -144};

double eps = 1e-3;

void Gaub_Jordan() {
	double matrix[N][N + 1];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N + 1; j++)
			matrix[i][j] = arr[i][j];

	
		double a;
		for (int k = 0; k < N; k++) {
			a = matrix[k][k];
			for (int j = k; j < N + 1; j++) {
				matrix[k][j] /= a;
			}

			for (int i = 0; i < N; i++) {
				if (i != k) {
					a = matrix[i][k];
					for (int j = k; j < N + 1; j++) {
						matrix[i][j] -= a * matrix[k][j];
					}
				}
			}
		}
		for (int i = 0; i < N; i++) {
			printf("%15.10lf", matrix[i][N]);
		}
		putchar('\n');
	
}

bool Direct_iter(double eps) {
	double matrix[N][N + 1];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N + 1; j++)
			matrix[i][j] = modern_arr[i][j];


	double x[N], oldx[N], currow, q = 0;
	for (int i = 0; i < N; i++) {
		currow = 0;
		for (int j = 0; j < N; j++) {
			if (i == j) continue;
			matrix[i][j] /= -matrix[i][i];
			currow += fabs(matrix[i][j]);
		}
		if (currow > q) q = currow;
		matrix[i][N] /= matrix[i][i];
		matrix[i][i] = 0;
		oldx[i] = matrix[i][N];
	}

	if (q >= 1) return false;

	bool flag = 1;
	while (flag) {
		flag = 0;
		for (int i = 0; i < N; i++) {
			x[i] = matrix[i][N];
			for (int j = 0; j < N; j++) {
				x[i] += matrix[i][j] * oldx[j];
			}
		}

		for (int i = 0; i < N; i++)
			if (fabs(oldx[i]) - fabs(x[i]) > (1 - q) / q * eps) flag = 1;
		for (int i = 0; i < N; i++)
			oldx[i] = x[i];
	}

	for (int i = 0; i < N; i++) {
		printf("%15.10lf", x[i]);
	}
	putchar('\n');

	return true;
}

int main() {

	for (int i = 0; i < N; i++)
		if (arr[i][i] == 0) {
			printf("matrix is not correct.\n");
			exit(1);
		}

	printf("Gaub_Jordan method:\n");
	Gaub_Jordan();
	printf("iteration method:\n");
	if (!Direct_iter(eps))
		 printf("Matrix is not for iteration method.\n");
	system("pause");
	return 0;
}
