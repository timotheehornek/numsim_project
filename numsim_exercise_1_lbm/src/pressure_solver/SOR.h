#pragma once

#include "pressure_solver.h"

#include <cmath>

class SOR
	: public Pressure_solver
{
private:
	const double m_w;
public:
	SOR(double dx, double dy, double eps, double max_it, double w);

	// run one  Successive OverRelaxation iteration and update p
	void run_it_step(Array2D& p, const Array2D& RHS, const std::array<int, 2> size) const override;
};