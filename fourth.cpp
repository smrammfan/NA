#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double func(double x) {
	return sqrt(3*x + 2);
}

double primitive_func(double x) {
	return (2./9)*sqrt(pow(3 * x + 2, 3));
}

double Integral(double n, double a, double b) {
	int i;
	double h; 
	double sig1 = 0.0; 
	double sig2 = 0.0; 

	n = 2 * ceil(n / 2);
	h = (b - a) / n;
	for (i = 1; i < n; i++) {
		if (i % 2 == 0)
			sig2 += func(a + i * h);
		else
			sig1 += func(a + i * h);
	}
	return h / 3 * (func(a) + func(b) + 4 * sig1 + 2 * sig2);
}

int main() {
	double n;
	double r, delta, h, maxY4, Io, In, I2n;
	double a = 1;
	double b = 17;
	double eps = 1e-4;

	
	Io = primitive_func(b) - primitive_func(a);
	printf("a =%.1f\n", a);
	printf("b =%.1f\n", b);
	printf("I = F(b) - F(a) = %f\n\n", Io);

	
	printf("Task1 - Simpson_method\n");
	printf("==================================================\n");
	printf("|    eps   |     h    |      I     |    delta    |\n");
	printf("==================================================\n");

	maxY4 = 8; 

	h = pow((180 * eps / (b - a) / maxY4), 0.25);
	n = 2 * ceil(0.5*(b - a) / h);
	h = (b - a) / n;
	delta = fabs(Io - Integral(n, a, b));
	printf("| %f | %.6f | %.8f | %.9f |\n", eps, h, Integral(n, a, b), delta);
	printf("==================================================\n\n");

	
	printf("Task2 Runge_mathod\n");
	printf("===========================================\n");
	printf("|    delta    |     h    |       Abs      |\n");
	printf("===========================================\n");
	n = ceil((b - a) / sqrt(sqrt(delta)));
	h = (b - a) / n;

	In = Integral(n, a, b);
	I2n = Integral(2 * n, a, b);
	r = fabs(In - I2n) / 15;

	while (r > delta) {
		n *= 2;
		In = I2n;
		I2n = Integral(2 * n, a, b);
		r = fabs(In - I2n) / 15;
	}

	printf("| %.9f | %.6f | %.12f |\n", delta, h, fabs(Io - I2n));
	printf("===========================================\n");
	system("pause");
	return 0;
}
