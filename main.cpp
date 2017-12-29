#include <iomanip>
#include <complex>
#include <cmath>
#include "Burgers.h"


int main()
{
	const std::complex<double> a_11 = (0.5, 0.5);
	int N = 200;
	int M = 300;
	double x_l = 0.0, x_r = 1.0;
	double T = 0.8;
	double eps = 0.1;

	//Создаем сетку
	Grid grd(N, M, x_l, x_r, T);
	double *x = grd.get_x_grid();
	double *t = grd.get_t_grid();

	//Создаем вектора граничных условий 
	double* u_l = new double[M + 1];
	double* u_r = new double[M + 1];
	for (int j = 0; j < (M + 1); j++) {
		*(u_l + j) = -5.0*cos(t[j]) + 8.0*sin(t[j]);
		*(u_r + j) = 2.0*cos(t[j]);
	}
	//Задаем начальную функцию u_0 
	double* u0 = new double[N + 1];
	for (int i = 0; i < (N + 1); i++)
		*(u0 + i) = ((x[i] + 1.0) + (x[i] - 5.0)*exp(-3.0*(x[i] - 0.50) / eps)) / (1.0 + exp(-3.0*(x[i] - 0.50) / eps));


	//Конструируем объект класса Burgers
	Burgers BurgersSol(N, M, x_l, x_r, T, a_11, u_l, u_r, u0, eps);
	delete[] u_l, u_r, u0;

	double** u = BurgersSol.get_solution();

	//Записываем в файл
	std::ofstream myofs;
	myofs.open("Results.txt");
	myofs << BurgersSol;

	for (int j = 0; j < (BurgersSol.get_M() + 1); j++)
		for (int i = 0; i < (BurgersSol.get_N() + 1); i++) {
			myofs << std::setprecision(4) << u[j][i] << std::endl;
		}

	myofs.close();

	return 0;
}