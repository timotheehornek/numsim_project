#include "SOR.h"

SOR::SOR(double eps, double max_it, double w)
	: Pressure_solver(eps, max_it), m_w{w} {}

void SOR::run_it_step(Discretization &discr,
					  const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
					  const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
					  const std::array<bool, 4> &useDirichletBc)
{
	for (int j{1}; j < discr.p().size()[1] - 1; ++j)
	{
		for (int i{1}; i < discr.p().size()[0] - 1; ++i)
		{
			if (!discr.is_in_obstacle(i, j, VAR_P))
				//! single node iteration step
				discr.p_ref(i, j) =
					(1 - m_w) * discr.p(i, j) + m_w																									   //< factor omega
													* std::pow(discr.dx() * discr.dy(), 2) / (2 * (std::pow(discr.dx(), 2) + std::pow(discr.dy(), 2))) //< prefactor
													* ((discr.p(i - 1, j) + discr.p(i + 1, j)) / std::pow(discr.dx(), 2)							   //< horizontal term
													   + (discr.p(i, j - 1) + discr.p(i, j + 1)) / std::pow(discr.dy(), 2)							   //< vertical term
													   - discr.RHS(i, j));																			   //< right hand side
		}
	}

	//! update boundary values around domain and obstacle
	update_boundaries(discr, bcBottom, bcTop, bcLeft, bcRight, useDirichletBc);
}

//! implementation of pressure solver
void SOR::solver(Discretization &discr,
							 const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
							 const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
							 const std::array<bool, 4> &useDirichletBc)
{
	//! get array size
	std::array<int, 2> size = discr.p().size();
	assert(discr.p().size() == discr.RHS().size());

	//! initialize iterations
	double res{residual(discr)};
	int it_counter{0};
	do
	{
		run_it_step(discr, bcBottom, bcTop, bcLeft, bcRight, useDirichletBc);
		res = residual(discr);
		++it_counter;
	} while (res > m_eps && it_counter < m_max_it);

	if (it_counter == m_max_it)
		std::cout << "max_it reached\n";
	DEBUG_PRINT(
		"Pressure solver terminated:\n"
		<< "\tCurrent iteration: " << it_counter << "\t\tCurrent residual: " << res << '\n'
		<< "\t   Max iterations: " << m_max_it << "\t    Max residual: " << m_eps);
}