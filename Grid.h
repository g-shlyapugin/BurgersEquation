#pragma once
class Grid {
public:
	Grid(int, int, double, double, double);
	~Grid();

	int get_N()const { return N; }
	int get_M()const { return M; }
	double get_T()const { return T; }

	double get_x_left()const { return x_left; }
	double get_x_right()const { return x_right; }

	double* get_x_grid()const { return x_grid; }
	double* get_t_grid()const { return t_grid; }



protected:
	int N, M;
	double x_left, x_right;
	double T;

	double* x_grid;
	double* t_grid;
};