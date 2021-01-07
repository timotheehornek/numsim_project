#include "pressure_solver.h"

#include <cmath>

Pressure_solver::Pressure_solver(double dx, double dy, double eps, double max_it, Message_passer_2D &MP)
	: m_dx{dx}, m_dy{dy}, m_eps{eps}, m_max_it{max_it}, m_MP{MP} 
	{}

double Pressure_solver::residual(const Staggered_grid &p, const Staggered_grid &RHS) const
{
	//! get array size
	std::array<int, 2> size = p.size();
	assert(p.size() == RHS.size());

	double p_xx{};
	double p_yy{};
	double residual{0.0};

	//! iterate over inner grid cells
	for (int j{1}; j < size[1] - 1; ++j)
	{
		for (int i{1}; i < size[0] - 1; ++i)
		{
			//! compute local result of laplace operator
			p_xx = (p(i - 1, j) - 2 * p(i, j) + p(i + 1, j)) / std::pow(m_dx, 2);
			p_yy = (p(i, j - 1) - 2 * p(i, j) + p(i, j + 1)) / std::pow(m_dy, 2);

			//! add squared local residual to result
			residual += std::pow(RHS(i, j) - (p_xx + p_yy), 2);
		}
	}
	//! compute number of inner cells
	const int N{(size[0] - 2) * (size[1] - 2)};
	//! return square root of sum of squared local residuals (L2 norm)
	//return std::pow(residual / N, .5);
	return residual / N;
}

//! implementation of pressure solver
void Pressure_solver::solver(Staggered_grid &p, const Staggered_grid &RHS)
{
	//! get array size
	std::array<int, 2> size = p.size();
	assert(p.size() == RHS.size());

	//! decide whether to run black-red or red-black
	bool black_red;
	int size_x_even{};
	int size_y_even{};
	std::array<int, 2> nCells_0_mod_2{};
	if (m_MP.rank() == 0)
	{
		size_x_even = size[0] % 2;
		size_y_even = size[1] % 2;
		//std::cout << "size_x_even= "<<size_x_even<<" size_y_even= "<<size_y_even<<'\n';
	}
	m_MP.broadcast(&size_x_even);
	m_MP.broadcast(&size_y_even);

	int p_x = m_MP.rank() % m_MP.prcs(0);
	int p_y = m_MP.rank() / m_MP.prcs(0);
	black_red = (size_x_even * p_x + size_y_even * p_y) % 2;

	//! initialize iterations
	double local_res{residual(p, RHS)};
	std::vector<double> all_local_res;
	all_local_res.resize(m_MP.size());

	double global_res{std::numeric_limits<double>::max()};

	int it_counter{0};

	while (global_res > m_eps && it_counter < m_max_it)
	{
		//! run local iteration step
		run_it_step(p, RHS, black_red);
		//! update local residual
		local_res = residual(p, RHS);
		//! update iteration counter
		++it_counter;
		m_MP.gather(local_res, all_local_res);
		if (m_MP.rank() == 0)
		{
			//! compute new global residual
			global_res = 0.0;
			for (double &n : all_local_res)
				global_res += n;
			global_res = std::pow(global_res / m_MP.size(), .5);
		}
		m_MP.broadcast(&global_res);
	}

	if (m_MP.rank() == 0)
	{
		DEBUG_PRINT(
			"Pressure solver terminated:\n"
			<< "\tCurrent iteration: " << it_counter << "\t\tCurrent residual: " << global_res << '\n'
			<< "\t   Max iterations: " << m_max_it << "\t    Max residual: " << m_eps);
	}
}