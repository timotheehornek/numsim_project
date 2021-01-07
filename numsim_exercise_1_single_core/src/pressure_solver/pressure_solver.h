#pragma once

#ifndef NDEBUG
#define DEBUG_PRINT(x) do { std::cout << x << '\n'; } while (0)
#else
#define DEBUG_PRINT(x) 
#endif

#include "array2d/array2d.h"

class Pressure_solver
{
protected:
	const double m_dx;
	const double m_dy;
	const double m_eps;
	const double m_max_it;

public:
	Pressure_solver(double dx, double dy, double eps, double max_it);
	
	//! run solver
	void solver(Array2D& p, const Array2D& RHS)const;

	//! run one iteration step
	virtual void run_it_step(Array2D& p, const Array2D& RHS, const std::array<int, 2> size) const = 0;

	// compute residual 
	double residual(const Array2D& p, const Array2D& RHS) const;
};