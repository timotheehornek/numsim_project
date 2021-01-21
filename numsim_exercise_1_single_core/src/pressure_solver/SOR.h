#pragma once

#include "pressure_solver.h"

#include <cmath>

class SOR
	: public Pressure_solver
{
private:
	const double m_w;

	enum SIDES
	{
		BOTTOM,
		TOP,
		LEFT,
		RIGHT
	};
	
public:
	SOR(double eps, double max_it, double w);

	// run one  Successive OverRelaxation iteration and update p
	//void run_it_step(Array2D& p, const Array2D& RHS, const std::array<int, 2> size) const override;
	void run_it_step(Discretization& discr,
	const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
	const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
	const std::array<bool, 4> &useDirichletBc) const override;
};