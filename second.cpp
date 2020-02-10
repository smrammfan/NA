#include <iostream>

using namespace std;
double function(double);
double first_derivative(double);
double second_derivative(double);
int roots_selection(double, double, double, double*);
double* method_iteration(double, double, double);
double*  method_hord(double, double, double);


int main() {
	double mas[4];
	double *res1 = (double*)malloc(sizeof(double) * 3);
	double *res2 = (double*)malloc(sizeof(double) * 3);
	int i, k;
	double eps, step = 0.3;
	double a, b;

	k = roots_selection(-10, 1, step, mas);
	
	printf("ITERATION:\n\n");
	for (i = 0; i <= k; i++) {
		a = mas[i];
		b = a + step;
		printf(" EPS            ROOT_%d              ACCURACY\n\n", i + 1);
		for (eps = 1e-2; eps >= 1e-14; eps *= 1e-3) {
			res1 = method_iteration(a, b, eps);
			printf("%0.E  %20.15f  %20.15f\n", eps, res1[0], res1[1]);
		}
		printf("\n");
	}

	printf("-----------------------------------------------------------\n");
	printf("HORD:\n\n");
	for (i = 0; i <= k; i++) {
		a = mas[i];
		b = a + step;
		printf(" EPS            ROOT_%d              ACCURACY\n\n", i + 1);
		for (eps = 1e-2; eps >= 1e-14; eps *= 1e-3) {
			res2 = method_hord(a, b, eps);
			printf("%0.E  %20.15f  %20.15f\n", eps, res2[0], res2[1]);
		}
		printf("\n");
	}

	printf("-----------------------------------------------------------\n");
	printf(" EPS      ITERATION      HORD \n\n");
	a = mas[0];
	b = a + step;
	for (eps = 1e-2; eps >= 1e-14; eps *= 1e-3) {
		res1 = method_iteration(a, b, eps);
		res2 = method_hord(a, b, eps);
		printf("%0.E   %6.0d     %6.0d\n", eps, (int)res1[2], (int)res2[2]);
	}
	
	system("PAUSE");


	return 0;
}

double function(double x){
	return exp(x)*sqrt(1 - x) - 0.1;
}

double first_derivative(double x) {
	return exp(x)*sqrt(1. - x) + exp(x)/(2*sqrt(1 - x));
}

double second_derivative(double x) {
	return  exp(x)*sqrt(1. - x) + exp(x) / (2 * sqrt(1 - x)) + ((exp(x)*(sqrt(1 - x) + 1)/(2*(1 - x))));

}

int roots_selection(double border1, double border2, double step, double *mas) {

	int i = -1;
	double a;

	for (a = border1; a <= border2; a += step)
		if ( function(a) * function(a + step) < 0)
			mas[++i] = a;
	return i;
}

double* method_iteration(double a, double b, double eps) {
	double m, M;
	double alfa, q;
	double x1, x;
	int n;
	double *res = (double*)malloc(sizeof(double) * 3);

	if (fabs(first_derivative(a)) > fabs(first_derivative(b))) {
		M = first_derivative(a);
		m = first_derivative(b);
	}
	else {
		M = first_derivative(b);
		m = first_derivative(a);
	}
	alfa = 1 / M;
	q = 1 - fabs(m / M);

	x = (b + a) / 2;
	x1 = x;
	n = 0;
	do {
		n++;
		x = x1;
		x1 = x - alfa * function(x);
	} while (fabs(x1 - x) > (1 - q) / q * eps);
	res[0] = x1;
	res[1] = fabs(fabs(x1) - fabs(x))*q / (1 - q);
	res[2] = (double)n;


	return res;
}

double*  method_hord(double a, double b, double eps) {
	double c, m, x, x1;
	int n;
	double *res = (double*)malloc(sizeof(double) * 3);

	m = fabs(first_derivative(a));
	for (x1 = a; x1 < b; x1 += 0.001) {
		if (fabs(first_derivative(x1)) < m)
			m = fabs(first_derivative(x1));
	}

	if (function(a) * second_derivative(a) > 0) {
		c = a;
		x1 = b;
	}
	else {
		c = b;
		x1 = a;
	}

	x = x1;
	n = 0;
	while (fabs(function(x)) / m > eps) {
		n++;
		x -= function(x)*(x - c) / (function(x) - function(c));
	}

	res[0] = x;
	res[1] = fabs(function(x)) / m;
	res[2] = (double)n;

	return res;
}
