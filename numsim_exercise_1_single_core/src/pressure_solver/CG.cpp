#include "CG.h"

CG::CG(double eps, int max_it, std::array<int, 2> nCells)
	: Pressure_solver(eps, max_it),
	  m_size{nCells[0] + 2, nCells[1] + 2},
	  m_residuals{nCells[0] + 2, nCells[1] + 2},
	  m_directions{nCells[0] + 2, nCells[1] + 2},
	  m_Ap{nCells[0] + 2, nCells[1] + 2}
{
}

//Matrix-vector-product
void CG::compute_Ap(Discretization &discr)
{
	// iterate over inner grid cells
	for (int j{1}; j < m_size[1] - 1; ++j)
	{
		for (int i{1}; i < m_size[0] - 1; ++i)
		{
			if (!discr.is_in_obstacle(i, j, VAR_P))
			{
				m_Ap(i, j) =
					((-2 * (std::pow(discr.dx(), 2) + std::pow(discr.dy(), 2))) / std::pow(discr.dx() * discr.dy(), 2)) * m_directions(i, j) //< central term
					+ (m_directions(i - 1, j) + m_directions(i + 1, j)) / std::pow(discr.dx(), 2)											 //< horizontal term
					+ (m_directions(i, j - 1) + m_directions(i, j + 1)) / std::pow(discr.dy(), 2);											 //< vertical term
			}
			else
			{
				m_Ap(i, j) = 0;
			}
		}
	}
}

// squared L2-norm of vector
double CG::compute_rk_transp_rk(const Staggered_grid &residuals)
{
	double norm_squared = 0;
	for (int j{0}; j < m_size[1]; ++j)
	{
		for (int i{0}; i < m_size[0]; ++i)
		{
			norm_squared += std::pow(residuals(i, j), 2); //! add local product to norm
		}
	}
	return norm_squared;
}

// squared L2-norm of p^T*A*p
double CG::compute_pk_transp_Apk()
{
	double result = 0;
	for (int j{0}; j < m_size[1]; ++j)
	{
		for (int i{0}; i < m_size[0]; ++i)
		{
			result += m_directions(i, j) * m_Ap(i, j); //! add local product to norm
		}
	}
	return result;
}

//! compute residual using p and RHS for each gridpoint; store value in "residuals"
void CG::compute_residual(Discretization &discr)
{
	double p_xx{};
	double p_yy{};

	//! iterate over inner grid cells
	for (int j{1}; j < m_size[1] - 1; ++j)
	{
		for (int i{1}; i < m_size[0] - 1; ++i)
		{
			if (!discr.is_in_obstacle(i, j, VAR_P))
			{
				//! compute local result of laplace operator
				p_xx = (discr.p(i - 1, j) - 2 * discr.p(i, j) + discr.p(i + 1, j)) / std::pow(discr.dx(), 2);
				p_yy = (discr.p(i, j - 1) - 2 * discr.p(i, j) + discr.p(i, j + 1)) / std::pow(discr.dy(), 2);

				//! add squared local residual to result
				m_residuals(i, j) = discr.RHS(i, j) - (p_xx + p_yy);
			}
			else
			{
				m_residuals(i, j) = 0;
			}
		}
	}
}

// compute new residual
void CG::update_residual(double alpha_k)
{

	for (int j{1}; j < m_size[1] - 1; ++j)
	{
		for (int i{1}; i < m_size[0] - 1; ++i)
		{
			// update local residual using auxiliary variables
			m_residuals(i, j) = m_residuals(i, j) - alpha_k * m_Ap(i, j);
		}
	}
}

// set new pressure values
void CG::update_p(double alpha_k, Discretization &discr)
{
	for (int j{1}; j < m_size[1] - 1; ++j)
	{
		for (int i{1}; i < m_size[0] - 1; ++i)
		{
			//! only update p for grid cells not in obstacle
			if (!discr.is_in_obstacle(i, j, VAR_P))
			{
				discr.p_ref(i, j) = discr.p_ref(i, j) + alpha_k * m_directions(i, j);
			}
		}
	}
}

// compute new search directions
void CG::update_directions(double beta_k)
{
	for (int j{1}; j < m_size[1] - 1; ++j)
	{
		for (int i{1}; i < m_size[0] - 1; ++i)
		{
			// update search directions using auxiliary variables such that they are conjugated to all previous search directions
			m_directions(i, j) = m_residuals(i, j) + beta_k * m_directions(i, j);
		}
	}
}

/*
//store all values of send in recv
void CG::copy_array(const Staggered_grid &send, Staggered_grid &recv, Discretization &discr)
{
	for (int j{1}; j < discr.p().size()[1] - 1; ++j)
	{
		for (int i{1}; i < discr.p().size()[0] - 1; ++i)
		{
			recv(i, j) = send(i, j);
		}
	}
}
*/
// run one  conjugated gradient iteration and update p
void CG::run_it_step(Discretization &discr)
{
	//! Using the most commonly used version of the CG-Algorithm
	//! Auxiliary variables (directions, residuals, Ap) are provided by Pressure_solver and already initialized

	compute_Ap(discr);
	double alpha_k = compute_rk_transp_rk(m_residuals) / compute_pk_transp_Apk();
	update_p(alpha_k, discr);
	Staggered_grid residual_k{m_size};
	//copy_array(residuals, residual_k, discr); //save old residuals
	residual_k = m_residuals;
	update_residual(alpha_k);
	double beta_k = compute_rk_transp_rk(m_residuals) / compute_rk_transp_rk(residual_k);
	update_directions(beta_k);
}

//! implementation of CG pressure solver
void CG::solver(Discretization &discr,
				const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
				const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
				const std::array<bool, 4> &useDirichletBc)
{
	//! get array size
	assert(discr.p().size() == discr.RHS().size() && discr.p().size() == m_size);

	//! initialize iterations
	double res{residual(discr)};
	int it_counter{0};

	//! initialize auxiliary variables
	compute_residual(discr);
	m_directions = m_residuals;

	//iterate until convergence (or max_it reached)
	do
	{
		run_it_step(discr);
		//res = residual(discr, residuals); //residuals are computed for each gridpoint in CG-algorithm, so build the global residual just with them
		res = residual(discr);
		++it_counter;
	} while (res > m_eps && it_counter < m_max_it);
	update_boundaries(discr, bcBottom, bcTop, bcLeft, bcRight, useDirichletBc); // set boundary conditions once AFTER iteration for each timestep is finished

	if (it_counter == m_max_it)
		std::cout << "max_it reached\n";
	DEBUG_PRINT(
		"Pressure solver terminated:\n"
		<< "\tCurrent iteration: " << it_counter << "\t\tCurrent residual: " << res << '\n'
		<< "\t   Max iterations: " << m_max_it << "\t    Max residual: " << m_eps);
}