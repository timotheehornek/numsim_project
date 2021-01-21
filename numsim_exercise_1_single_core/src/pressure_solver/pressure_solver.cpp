#include "pressure_solver.h"


Pressure_solver::Pressure_solver(double eps, double max_it, const std::array<bool,4>& p_0_boundary)
	: m_eps{ eps }, m_max_it{ max_it }, m_p_0_boundary{p_0_boundary}{}


//double Pressure_solver::residual(const Array2D& p, const Array2D& RHS) const
	double Pressure_solver::residual(Discretization& discr) const
{
	//! get array size
	std::array<int, 2>size = discr.p().size();
	assert(discr.p().size() == discr.RHS().size());

	double p_xx{};
	double p_yy{};
	double residual{ 0.0 };

	//! iterate over inner grid cells
	for (int j{ 1 }; j < size[1] - 1; ++j)
	{
		for (int i{ 1 }; i < size[0] - 1; ++i)
		{
			if (!discr.is_in_obstacle(i,j,VAR_P))
			{
				//! compute local result of laplace operator
				p_xx = (discr.p(i - 1, j) - 2 * discr.p(i, j) + discr.p(i + 1, j)) / std::pow(discr.dx(), 2);
				p_yy = (discr.p(i, j - 1) - 2 * discr.p(i, j) + discr.p(i, j + 1)) / std::pow(discr.dy(), 2);
				
				//! add squared local residual to result
				residual += std::pow(discr.RHS(i,j) - (p_xx + p_yy), 2);
			}
		}
	}
	//! compute number of inner cells
	const int N{ (size[0] - 2) * (size[1] - 2) };
	//! return square root of sum of squared local residuals (L2 norm)
	return std::pow(residual / N, .5);
}

//! implementation of pressure solver
//void Pressure_solver::solver(Array2D& p, const Array2D& RHS) const
void Pressure_solver::solver(Discretization& discr) const
{
	//! get array size
	std::array<int, 2>size = discr.p().size();
	assert(discr.p().size() == discr.RHS().size());

	//! initialize iterations
	double res{ residual(discr) };
	int it_counter{ 0 };
	do
	{
		run_it_step(discr);
		res = residual(discr);
		++it_counter;
	} while (res > m_eps && it_counter < m_max_it);

	if (it_counter == m_max_it)
		std::cout << "max_it reached\n";
	DEBUG_PRINT
	(
		"Pressure solver terminated:\n"
		<< "\tCurrent iteration: " << it_counter << "\t\tCurrent residual: " << res << '\n'
		<< "\t   Max iterations: " << m_max_it << "\t    Max residual: " << m_eps
	);
}