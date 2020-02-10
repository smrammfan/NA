#include <iostream>
#include <cmath>

const double ln2 = 0.693147180559945309417232121458;
using namespace std;

void drawTable();
void calculate(double, double, int);

int main() {

	drawTable();
	system("pause");
	return 0;
}

void drawTable() {
	double startEPS = 1e-2;
	double endEPS = 1e-14;
	double step = 1e-3;
	double a = 0.1;
	double b = 15.0;
	double h = (b - a) / 10.0;
	double x = (a + b) / 2;

	//task 1;
	printf("---------------------------------------------------------------------\n");
	printf("|%15s |%5s |%20s |%20s |\n", "eps", "n", "delta", "R");
	printf("---------------------------------------------------------------------\n");
	for (double d = startEPS; d >= endEPS; d *= step) {
		calculate(d, x, 1);
	}
	printf("\n\n");
	//task 2
	printf("---------------------------------------------------------\n");
	printf("|%10s |%20s |%20s |\n", "X", "delta", "R");
	printf("---------------------------------------------------------\n");

	x = a;
	for (int i = 1; i <= 11; i++) {
		calculate(1e-8, x, 2);
		x = a + h * i;
	}
}

void calculate(double eps, double x, int task) {

	int m = 0;

	double z = x;

	for (; z >= 1; m++) {
		z /= 2;
	}

	for (; z <= 0.5; m--) {
		z *= 2;
	}

	double a = (1 - z) / (1 + z);
	int n = 1;
	double L = 1;
	double sumL = 0;
	double fourEPS = 4 * eps;

	while (!(L < fourEPS)) {
		L = pow(a, 2 * n - 1) / (2 * n - 1);
		n++;
		sumL += L;
		
	}

	double R = L / 4;
	double ln = m * ln2 - 2 * sumL - R;
	double delta = abs(ln - log(x));

	if (task == 1)
		printf("|%15E |%5d |%20E |%20E |\n", eps, n, delta, R);
	if (task == 2)
		printf("|%10f |%20E |%20E |\n", x, delta, R);
}

