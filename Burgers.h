#pragma once
#include "Grid.h"
#include <fstream> 
class Burgers : public Grid {
private:
	std::complex<double> a_11;		//Параметр схемы CROS1а
	double* u_left;					//Левое граничное условие
	double* u_right;				//Правое граничное условие
	double* u_0;					//Начальная функция
	double eps;						//Малый (сингулярный) параметр 
	double** u_sol;
public:
	Burgers(int n, int m, double x_l, double x_r,
		double T1, std::complex<double> a, double* ul,
		double*, double* u_00, double ep);

	//В случае, если граничные условия не зависят от времени
	Burgers(int n, int m, double x_l, double x_r,
		double T1, std::complex<double> a,
		double ul, double ur, double* u_00, double ep);

	double** get_solution();

	//Деструктор
	~Burgers();

};

std::ofstream& operator<<(std::ofstream& ofs, const Burgers& burg_obj);

//Функция, реализующая метод прогонки
std::complex <double>* TridiagonalMatrixAlgorithm(int N, std::complex<double>* a,
	std::complex<double>* b, std::complex<double>* c, double* d);

//Функ, вычисляющая правие части системы ОДУ
double* fCalculation(int N, double* u, double h, double eps, double u_left, double u_right);

//Вычисление диагоналей 3-х диаг. матр. на каждом шаге
void DiagonalsPreparation(std::complex <double>* a, std::complex<double>* b, std::complex<double>* c,
	int N, double* u, double tau, double h, double eps, double u_left, double u_right, std::complex<double> a_11);
