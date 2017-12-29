#include "Grid.h"

Grid::Grid(int n, int m, double x_l, double x_r, double t) :
	N(n), M(m), x_left(x_l), x_right(x_r), T(t) {
	// Шаг сетки по координате и времени
	double h = (x_right - x_left) / N;
	double tau = T / M;

	//Формируем сетки по координате и времени												
	x_grid = new double[N + 1];
	for (int i = 0; i < (N + 1); i++)
		x_grid[i] = x_left + h * i;

	t_grid = new double[M + 1];
	for (int i = 0; i < (M + 1); i++)
		t_grid[i] = 0.0 + tau * i;
}

Grid::~Grid() {
	delete[] x_grid;
	delete[] t_grid;
}