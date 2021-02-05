#include "pressure_solver.h"

Pressure_solver::Pressure_solver(double eps, double max_it)
	: m_eps{eps}, m_max_it{max_it} {}

double Pressure_solver::residual(Discretization &discr) const
{
	//! get array size
	std::array<int, 2> size = discr.p().size();
	assert(discr.p().size() == discr.RHS().size());

	double p_xx{};
	double p_yy{};
	double residual{0.0};

	//! iterate over inner grid cells
	for (int j{1}; j < size[1] - 1; ++j)
	{
		for (int i{1}; i < size[0] - 1; ++i)
		{
			if (!discr.is_in_obstacle(i, j, VAR_P))
			{
				//! compute local result of laplace operator
				p_xx = (discr.p(i - 1, j) - 2 * discr.p(i, j) + discr.p(i + 1, j)) / std::pow(discr.dx(), 2);
				p_yy = (discr.p(i, j - 1) - 2 * discr.p(i, j) + discr.p(i, j + 1)) / std::pow(discr.dy(), 2);

				//! add squared local residual to result
				residual += std::pow(discr.RHS(i, j) - (p_xx + p_yy), 2);
			}
		}
	}
	//! compute number of inner cells
	const int N{(size[0] - 2) * (size[1] - 2)};
	//! return square root of sum of squared local residuals (L2 norm)
	return std::pow(residual / N, .5);
}

void Pressure_solver::update_boundaries(Discretization& discr,
    const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
    const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
    const std::array<bool, 4> &useDirichletBc) const
{
  //! update boundary values on left and right
	if (useDirichletBc[LEFT]) //< Neuman for pressure
		for (int j{1}; j < discr.p().size()[1] - 1; ++j)
			discr.p_ref(0, j) = discr.p(1, j);
	else //< Dirichlet for pressure
		for (int j{1}; j < discr.p().size()[1] - 1; ++j)
			discr.p_ref(0, j) = 2 * bcLeft[0] - discr.p(1, j);
	if (useDirichletBc[RIGHT]) //< Neuman for pressure
		for (int j{1}; j < discr.p().size()[1] - 1; ++j)
			discr.p_ref(discr.p().size()[0] - 1, j) = discr.p(discr.p().size()[0] - 2, j);
	else //< Dirichlet for pressure
		for (int j{1}; j < discr.p().size()[1] - 1; ++j)
			discr.p_ref(discr.p().size()[0] - 1, j) = 2 * bcRight[0] - discr.p(discr.p().size()[0] - 2, j);

	//! update boundary values on bottom and top
	if (useDirichletBc[BOTTOM]) //< Neuman for pressure
		for (int i{0}; i < discr.p().size()[0]; ++i)
			discr.p_ref(i, 0) = discr.p(i, 1);
	else //< Dirichlet for pressure
		for (int i{0}; i < discr.p().size()[0]; ++i)
			discr.p_ref(i, 0) = 2 * bcBottom[1] - discr.p(i, 1);
	if (useDirichletBc[TOP]) //< Neuman for pressure
		for (int i{0}; i < discr.p().size()[0]; ++i)
			discr.p_ref(i, discr.p().size()[1] - 1) = discr.p(i, discr.p().size()[1] - 2);
	else //< Dirichlet for pressure
		for (int i{0}; i < discr.p().size()[0]; ++i)
			discr.p_ref(i, discr.p().size()[1] - 1) = 2 * bcTop[1] - discr.p(i, discr.p().size()[1] - 2);

	//! update boundary values around obstacle
	if (discr.obstacle_exist())
	{
		// top + bottom (without corners)
		for (int i{discr.obstacle_pos(0) + 2}; i <= discr.obstacle_pos(2); ++i)
		{
			// top
			discr.p_ref(i, discr.obstacle_pos(3) + 1) = discr.p(i, discr.obstacle_pos(3) + 2);

			// bottom
			discr.p_ref(i, discr.obstacle_pos(1) + 1) = discr.p(i, discr.obstacle_pos(1));
		}

		// left + right (without corners)
		for (int j{discr.obstacle_pos(1) + 2}; j <= discr.obstacle_pos(3); ++j)
		{
			// left
			discr.p_ref(discr.obstacle_pos(0) + 1, j) = discr.p(discr.obstacle_pos(0), j);

			// right
			discr.p_ref(discr.obstacle_pos(2) + 1, j) = discr.p(discr.obstacle_pos(2) + 2, j);
		}

		// boundary corners
		discr.p_ref(discr.obstacle_pos(0) + 1, discr.obstacle_pos(1) + 1) = .5 * (discr.p(discr.obstacle_pos(0), discr.obstacle_pos(1) + 1) + discr.p(discr.obstacle_pos(0) + 1, discr.obstacle_pos(1)));
		discr.p_ref(discr.obstacle_pos(2) + 1, discr.obstacle_pos(1) + 1) = .5 * (discr.p(discr.obstacle_pos(2) + 2, discr.obstacle_pos(1) + 1) + discr.p(discr.obstacle_pos(2) + 1, discr.obstacle_pos(1)));
		discr.p_ref(discr.obstacle_pos(2) + 1, discr.obstacle_pos(3) + 1) = .5 * (discr.p(discr.obstacle_pos(2) + 2, discr.obstacle_pos(3) + 1) + discr.p(discr.obstacle_pos(2) + 1, discr.obstacle_pos(3) + 2));
		discr.p_ref(discr.obstacle_pos(0) + 1, discr.obstacle_pos(3) + 1) = .5 * (discr.p(discr.obstacle_pos(0), discr.obstacle_pos(3) + 1) + discr.p(discr.obstacle_pos(0) + 1, discr.obstacle_pos(3) + 2));
	}
}