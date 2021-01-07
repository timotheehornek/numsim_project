#pragma once

#include "pressure_solver.h"

#include <cmath>

class Gauss_seidel :
	public Pressure_solver
{
public:
	Gauss_seidel(double dx, double dy, double eps, double max_it);

	// run one Gauss-Seidel iteration and update p
	void run_it_step(Array2D& p, const Array2D& RHS, const std::array<int, 2> size) const override;
};
