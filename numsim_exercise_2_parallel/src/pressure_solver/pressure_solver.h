#pragma once

#ifndef NDEBUG
#define DEBUG_PRINT(x) do { std::cout << x << '\n'; } while (0)
#else
#define DEBUG_PRINT(x) 
#endif

#include "array2d/array2d.h"
#include "discretization/discretization.h"
#include "discretization/staggered_grid.h"

#include <algorithm>
#include <limits>
#include <cmath>
#include <vector>

class Pressure_solver
{
protected:
	const double m_dx;
	const double m_dy;
	const double m_eps;
	const double m_max_it;

	Message_passer_2D &m_MP;
	

public:
	Pressure_solver(double dx, double dy, double eps, double max_it, Message_passer_2D &MP);
	
	//! run solver
	//void solver(Array2D& p, const Array2D& RHS)const;
	void solver(Staggered_grid& p, const Staggered_grid& RHS);

	//! run one iteration step
	virtual void run_it_step(Staggered_grid &p, const Staggered_grid &RHS, bool black_red) = 0;

	// compute residual 
	double residual(const Staggered_grid& p, const Staggered_grid& RHS) const;
};