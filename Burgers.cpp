#include <complex>
#include <cmath>
#include "Burgers.h"

Burgers::Burgers(int n, int m, double x_l, double x_r,
	double T1, std::complex<double> a, double* ul, double* ur,
	double* u_00, double ep) :
	Grid(n, m, x_l, x_r, T1), a_11(a), eps(ep) {
	//
	u_left = new double[m + 1];
	u_right = new double[m + 1];
	for (int i = 0; i < (m + 1); i++) {
		*(u_left + i) = *(ul + i);
		*(u_right + i) = *(ur + i);
	}

	u_0 = new double[n + 1];
	for (int i = 0; i < (n + 1); i++)
		*(u_0 + i) = *(u_00 + i);

	//Формируем массив решений
	u_sol = new double*[M + 1];
	for (int count = 0; count < (M + 1); count++)
		u_sol[count] = new double[N + 1];
}


//Постоянные ГУ
Burgers::Burgers(int n, int m, double x_l, double x_r,
	double T1, std::complex<double> a, double ul, double ur, double* u_00, double ep) :
	Grid(n, m, x_l, x_r, T1), a_11(a), eps(ep) {
	//
	u_left = new double[m + 1];
	u_right = new double[m + 1];
	for (int i = 0; i < (m + 1); i++) {
		*(u_left + i) = ul;
		*(u_right + i) = ur;
	}

	u_0 = new double[n + 1];
	for (int i = 0; i < (n + 1); i++)
		*(u_0 + i) = *(u_00 + i);

	//Формируем массив решений
	u_sol = new double*[M + 1];
	for (int count = 0; count < (M + 1); count++)
		u_sol[count] = new double[N + 1];
}

double** Burgers::get_solution() {
	//Диагонали трехдиагональной матрицы
	std::complex<double> *a = new std::complex<double>[N - 1];
	std::complex<double> *b = new std::complex<double>[N - 1];
	std::complex<double> *c = new std::complex<double>[N - 1];
	for (int count = 0; count < (N - 1); count++) {
		a[count] = (0.0, 0.0);
		b[count] = (0.0, 0.0);
		c[count] = (0.0, 0.0);
	}

	//Вектор правых частей
	double *f = new double[N - 1];
	for (int count = 0; count < (N - 1); count++)
		f[count] = 0.0;

	//Вспомогательные вектор 
	std::complex<double> *w = new std::complex<double>[N - 1];
	for (int count = 0; count < (N - 1); count++)
		w[count] = (0.0, 0.0);
	//

	//Задаём начальное условие
	for (int i = 0; i < (N + 1); i++)
		u_sol[0][i] = u_0[i];

	//Шаги сетки
	double h = (x_right - x_left) / N;
	double tau = T / M;

	//Ищем решение ОДУ в каждый следующий момент времени (t + tau) по известному решению в текущй момент времени t
	for (int j = 0; j < M; j++) {
		//Cхема Розенброка(CROS1)
		//Подготавливаем массивы, которые будут содержать коэффициенты на диагоналях трёхдиагональной матрицы
		DiagonalsPreparation(a, b, c, N, u_sol[j], tau, h, eps, u_left[j], u_right[j], a_11);

		//Вычисляем вектор правой части решаемой системы ОДУ
		f = fCalculation(N, u_sol[j], h, eps, u_left[j], u_right[j]);

		//Решаем систему линейных уравнений методом прогонки
		w = TridiagonalMatrixAlgorithm(N, a, b, c, f);

		for (int i = 0; i < (N - 1); i++) {
			//Вычисляем действительную часть w
			double ww = w[i].real();
			//Решение на следующем временном слое
			u_sol[j + 1][i + 1] = u_sol[j][i + 1] + (t_grid[j + 1] - t_grid[j])*ww;
		}

		//Граничные условия
		u_sol[j + 1][0] = u_left[j + 1];
		u_sol[j + 1][N] = u_right[j + 1];
	}

	delete[] a; delete[] b; delete[] c; delete[] f; delete[] w;
	return u_sol;
}


Burgers::~Burgers() {
	delete[] u_0;
	delete[] u_right;
	delete[] u_left;
	for (int count = 0; count < (M + 1); count++)
		delete[] u_sol[count];
}

std::ofstream& operator<<(std::ofstream& ofs, const Burgers& burg_obj) {
	ofs << burg_obj.get_N() << std::endl << burg_obj.get_M() << std::endl << burg_obj.get_x_left() <<
		std::endl << burg_obj.get_x_right() << std::endl << burg_obj.get_T() << std::endl;
	return ofs;
}


std::complex<double>* TridiagonalMatrixAlgorithm(int N, std::complex <double>* a,
	std::complex<double>* b, std::complex<double>* c, double* d) {
	std::complex<double> *v = new std::complex<double>[N - 1];
	for (int count = 0; count < (N - 1); count++)
		v[count] = (0.0, 0.0);

	std::complex<double> *x = new std::complex<double>[N - 1];
	for (int count = 0; count < (N - 1); count++)
		x[count] = (0.0, 0.0);

	std::complex<double> w = a[0];
	x[0] = d[0] / w;
	for (int i = 1; i < (N - 1); i++) {
		v[i - 1] = c[i - 1] / w;
		w = a[i] - b[i] * v[i - 1];
		x[i] = (d[i] - b[i] * x[i - 1]) / w;
	}

	for (int i = (N - 3); i >= 0; i--)
		x[i] = x[i] - v[i] * x[i + 1];

	delete[] v;

	return x;
}

double* fCalculation(int N, double* u, double h, double eps, double u_left, double u_right) {
	//Функция, которая вычисляет вектор правой части решаемой системы ОДУ
	double *f = new double[N - 1];
	for (int count = 0; count < (N - 1); count++)
		f[count] = (0.0, 0.0);
	//Учитываем левое граничное условие отдельно от общего цикла
	f[0] = (eps / (h*h))*(u_left - 2 * u[1] + u[2]) - u[1] * (u[2] - u_left) / (2 * h);
	// Общий цикл
	for (int i = 1; i < (N - 2); i++)
		f[i] = (eps / (h*h))*(u[i] - 2 * u[i + 1] + u[i + 2]) - u[i + 1] * (u[i + 2] - u[i]) / (2 * h);
	//Учитываем правое граничное условие отдельно от общего цикла
	f[N - 2] = (eps / (h*h))*(u[N - 2] - 2 * u[N - 1] + u_right) - u[N - 1] * (u_right - u[N - 2]) / (2 * h);

	return f;
}

void DiagonalsPreparation(std::complex<double>* a, std::complex<double>* b, std::complex<double>* c,
	int N, double* u, double tau, double h, double eps, double u_left, double u_right, std::complex<double> a_11) {
	//Функция, которая подготавливает массивы, которые содержат диагонали матрицы решаемой системы ОДУ.
	//Данная матрица имеет вид[E - a_11*tau*f_u(u, t)] и является трёхдиагональной a -- главная диагональ, 
	// b --нижняя, c--верхняя

	a[0] = 1.0 - a_11 * tau*(-2 * eps / (h*h) - (u[2] - u_left) / (2 * h));
	c[0] = -a_11 * tau*(eps / (h*h) - u[1] / (2 * h));

	for (int i = 1; i < (N - 2); i++) {
		b[i] = -a_11 * tau*(eps / (h*h) + u[i + 1] / (2 * h));
		a[i] = 1.0 - a_11 * tau*(-2 * eps / (h*h) - (u[i + 2] - u[i]) / (2 * h));
		c[i] = -a_11 * tau*(eps / (h*h) - u[i + 1] / (2 * h));
	}

	b[N - 2] = -a_11 * tau*(eps / (h*h) + u[N - 1] / (2 * h));
	a[N - 2] = 1.0 - a_11 * tau*(-2 * eps / (h*h) - (u_right - u[N - 2]) / (2 * h));
}