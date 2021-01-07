#include "gauss_seidel.h"

Gauss_seidel::Gauss_seidel(double dx, double dy, double eps, double max_it)
	: Pressure_solver(dx,dy,eps,max_it)
{}

void Gauss_seidel::run_it_step(Array2D& p, const Array2D& RHS, const std::array<int, 2> size) const
{
	for (int j{ 1 }; j < size[1] - 1; ++j)
	{
		for (int i{ 1 }; i < size[0] - 1; ++i)
		{
			//! single node iteration step
			p(i, j) =
				std::pow(m_dx * m_dy, 2) / (2 * (std::pow(m_dx, 2) + std::pow(m_dy, 2))) //< prefactor
				* ((p(i - 1, j) + p(i + 1, j)) / std::pow(m_dx, 2)						 //< horizontal term
					+ (p(i, j - 1) + p(i, j + 1)) / std::pow(m_dy, 2)					 //< vertical term
					- RHS(i, j));														 //< right hand side
		}
		//! update boundary values on left and right
		p(0, j) = p(1, j);
		p(size[0] - 1, j) = p(size[0] - 2, j);
	}

	//! update boundary values on bottom and top
	for (int i{ 0 }; i < size[0] ; ++i)
	{
		p(i, 0) = p(i, 1);
		p(i, size[1] - 1) = p(i, size[1] - 2);
	}
}
