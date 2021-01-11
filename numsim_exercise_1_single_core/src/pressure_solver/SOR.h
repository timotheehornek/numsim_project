#pragma once

#include "pressure_solver.h"

#include <cmath>

class SOR
	: public Pressure_solver
{
private:
	const double m_w;
	
public:
	SOR(double eps, double max_it, double w, const std::array<bool,4>& p_0_boundary);

	// run one  Successive OverRelaxation iteration and update p
	//void run_it_step(Array2D& p, const Array2D& RHS, const std::array<int, 2> size) const override;
	void run_it_step(Discretization& discr) const override;
};