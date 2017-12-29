#pragma once
#include "Grid.h"
#include <fstream> 
class Burgers : public Grid {
private:
	std::complex<double> a_11;		//�������� ����� CROS1�
	double* u_left;					//����� ��������� �������
	double* u_right;				//������ ��������� �������
	double* u_0;					//��������� �������
	double eps;						//����� (�����������) �������� 
	double** u_sol;
public:
	Burgers(int n, int m, double x_l, double x_r,
		double T1, std::complex<double> a, double* ul,
		double*, double* u_00, double ep);

	//� ������, ���� ��������� ������� �� ������� �� �������
	Burgers(int n, int m, double x_l, double x_r,
		double T1, std::complex<double> a,
		double ul, double ur, double* u_00, double ep);

	double** get_solution();

	//����������
	~Burgers();

};

std::ofstream& operator<<(std::ofstream& ofs, const Burgers& burg_obj);

//�������, ����������� ����� ��������
std::complex <double>* TridiagonalMatrixAlgorithm(int N, std::complex<double>* a,
	std::complex<double>* b, std::complex<double>* c, double* d);

//����, ����������� ������ ����� ������� ���
double* fCalculation(int N, double* u, double h, double eps, double u_left, double u_right);

//���������� ���������� 3-� ����. ����. �� ������ ����
void DiagonalsPreparation(std::complex <double>* a, std::complex<double>* b, std::complex<double>* c,
	int N, double* u, double tau, double h, double eps, double u_left, double u_right, std::complex<double> a_11);
