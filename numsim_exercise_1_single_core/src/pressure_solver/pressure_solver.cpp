#include "pressure_solver.h"


Pressure_solver::Pressure_solver(double dx, double dy, double eps, double max_it)
	: m_dx{ dx }, m_dy{ dy }, m_eps{ eps }, m_max_it{ max_it }{}


double Pressure_solver::residual(const Array2D& p, const Array2D& RHS) const
{
	//! get array size
	std::array<int, 2>size = p.size();
	assert(p.size() == RHS.size());

	double p_xx{};
	double p_yy{};
	double residual{ 0.0 };

	//! iterate over inner grid cells
	for (int j{ 1 }; j < size[1] - 1; ++j)
	{
		for (int i{ 1 }; i < size[0] - 1; ++i)
		{
			//! compute local result of laplace operator
			p_xx = (p(i - 1, j) - 2 * p(i, j) + p(i + 1, j)) / std::pow(m_dx, 2);
			p_yy = (p(i, j - 1) - 2 * p(i, j) + p(i, j + 1)) / std::pow(m_dy, 2);

			//! add squared local residual to result
			residual += std::pow(RHS(i,j) - (p_xx + p_yy), 2);
		}
	}
	//! compute number of inner cells
	const int N{ (size[0] - 2) * (size[1] - 2) };
	//! return square root of sum of squared local residuals (L2 norm)
	return std::pow(residual / N, .5);
}

//! implementation of pressure solver
//void Pressure_solver::solver(Array2D& p, const Array2D& RHS) const
void Pressure_solver::solver(Discretization& discr) const;
{
	//! get array size
	std::array<int, 2>size = discr.p().size();
	assert(discr.p().size() == discr.RHS().size());

	//! initialize iterations
	double res{ residual(discr.p(),discr.RHS()) };
	int it_counter{ 0 };
	while (res > m_eps && it_counter < m_max_it)
	{
		run_it_step(discr);
		res = residual(discr.p(), discr.RHS());
		++it_counter;
	}

	DEBUG_PRINT
	(
		"Pressure solver terminated:\n"
		<< "\tCurrent iteration: " << it_counter << "\t\tCurrent residual: " << res << '\n'
		<< "\t   Max iterations: " << m_max_it << "\t    Max residual: " << m_eps
	);
}