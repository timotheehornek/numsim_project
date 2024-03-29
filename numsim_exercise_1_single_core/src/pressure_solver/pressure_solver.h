#pragma once

#ifndef NDEBUG
#define DEBUG_PRINT(x)          \
	do                          \
	{                           \
		std::cout << x << '\n'; \
	} while (0)
#else
#define DEBUG_PRINT(x)
#endif

#include "array2d/array2d.h"
#include "discretization/discretization.h"

#include <array>
#include <cmath>

class Pressure_solver
{
protected:
	enum VARS
	{
		VAR_U,
		VAR_V,
		VAR_P
	};
	enum SIDES
	{
		BOTTOM,
		TOP,
		LEFT,
		RIGHT
	};

	const double m_eps;
	const double m_max_it;

	//! compute residual
	double residual(Discretization &discr) const;

public:
	Pressure_solver(double eps, double max_it);

	//! run solver
	virtual void solver(Discretization &discr,
						const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
						const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
						const std::array<bool, 4> &useDirichletBc) = 0;
	//! update boundary values of p and store them
	void update_boundaries(Discretization &discr,
						   const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
						   const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
						   const std::array<bool, 4> &useDirichletBc) const;
};
