#include "discretization.h"

Discretization::Discretization(const std::array<int, 2> &nCells, const std::array<double, 2> &physicalSize, const double re, const std::array<double, 2> &g)
	: m_nCells{nCells},
	  m_u{nCells[0] + 1, nCells[1] + 2},
	  m_v{nCells[0] + 2, nCells[1] + 1},
	  m_p{nCells[0] + 2, nCells[1] + 2},
	  m_F{nCells[0] + 2, nCells[1] + 2},
	  m_G{nCells[0] + 2, nCells[1] + 2},
	  m_RHS{nCells[0] + 2, nCells[1] + 2},
	  m_dx{physicalSize[0] / nCells[0]},
	  m_dy{physicalSize[1] / nCells[1]},
	  m_re{re},
	  m_g{g}
{
}

//! compute dt from stability conditions
void Discretization::set_dt(double dt_max, double sec_factor)
{
	//! stability condition induced by diffusion operator
	const double stab_cond_1{
		m_re / 2 *
		(std::pow(m_dx, 2) * std::pow(m_dy, 2)) / (std::pow(m_dx, 2) + std::pow(m_dy, 2))};

	//! stability conditions induced by convection operator
	const double stab_cond_2{m_dx / std::abs(m_u.abs_max())};
	const double stab_cond_3{m_dy / std::abs(m_v.abs_max())};

	//! find largest possible dt considering security factor
	const double dt_best{std::min(std::min(stab_cond_1, stab_cond_2), stab_cond_3) * sec_factor};

	//! set dt respecting user input for highest dt
	m_dt = std::min(dt_best, dt_max);

	assert(m_dt > 0);
}

//! set dt by taking dt from user
void Discretization::set_dt(double dt_max)
{
	//! set dt respecting user input for highest dt
	m_dt = dt_max;

	assert(m_dt > 0);
}

//! getter
double Discretization::dt() const
{
	return m_dt;
}
double Discretization::dx() const
{
	return m_dx;
}
double Discretization::dy() const
{
	return m_dy;
}
const Staggered_grid &Discretization::u() const
{
	return m_u;
}
const Staggered_grid &Discretization::v() const
{
	return m_v;
}
const Staggered_grid &Discretization::p() const
{
	return m_p;
}
const Staggered_grid &Discretization::F() const
{
	return m_F;
}
const Staggered_grid &Discretization::G() const
{
	return m_G;
}
const Staggered_grid &Discretization::RHS() const
{
	return m_RHS;
}
Staggered_grid &Discretization::p_ref()
{
	return m_p;
}
const std::array<int, 2> &Discretization::nCells() const
{
	return m_nCells;
}

void Discretization::setup_bound_val_uv(
	const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
	const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
	const std::array<bool, 4> &useDirichletBc)
{
	enum SIDES
	{
		BOTTOM,
		TOP,
		LEFT,
		RIGHT
	};

	// bottom
	if (useDirichletBc[BOTTOM]) //< Dirichlet
	{
		for (int i{1}; i < m_u.x_max() - 1; ++i)
			m_u(i, 0) = 2 * bcBottom[0] - m_u(i, 1);
		for (int i{0}; i < m_v.x_max(); ++i)
			m_v(i, 0) = bcBottom[1];
	}
	else //< Neumann
	{
		//m_v(i, 0) = - bcBottom[1] * m_dy;
	}
	//####################
	//HIER WEITERMACHEN

	// top
	if (useDirichletBc[BOTTOM])

	//! set u
	for (int j{0}; j < m_u.y_max(); ++j)
	{
		//! left boundary
		m_u(0, j) = bcLeft[0];
		//! right boundary
		m_u(m_u.x_max() - 1, j) = bcRight[0];
	}
	for (int i{1}; i < m_u.x_max() - 1; ++i)
	{
		//! bottom boundary
		m_u(i, 0) = 2 * bcBottom[0] - m_u(i, 1);
		//! top boundary
		m_u(i, m_u.y_max() - 1) = 2 * bcTop[0] - m_u(i, m_u.y_max() - 2);
	}

	//! set v
	for (int j{1}; j < m_v.y_max() - 1; ++j)
	{
		//! left boundary
		m_v(0, j) = 2 * bcLeft[1] - m_v(1, j);
		//! right boundary
		m_v(m_v.x_max() - 1, j) = 2 * bcRight[1] - m_v(m_v.x_max() - 2, j);
	}
	for (int i{0}; i < m_v.x_max(); ++i)
	{
		// bottom boundary
		m_v(i, 0) = bcBottom[1];
		// top boundary
		m_v(i, m_v.y_max() - 1) = bcTop[1];
	}
}

void Discretization::update_bound_val_uv(
	const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
	const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
	const std::array<bool, 4> &useDirichletBc)
{
	enum SIDES
	{
		BOTTOM,
		TOP,
		LEFT,
		RIGHT
	};

	// bottom
	if (useDirichletBc[BOTTOM]) //< Dirichlet
		for (int i{1}; i < m_u.x_max() - 1; ++i)
			m_u(i, 0) = 2 * bcBottom[0] - m_u(i, 1);
	else //< Neumann
		for (int i{1}; i < m_u.x_max() - 1; ++i)
		{
			m_u(i, 0) = m_u(i, 1);
			m_v(i, 0) = m_v(i, 1) - bcBottom[1] * m_dy;
		}

	// top
	if (useDirichletBc[TOP]) //< Dirichlet
		for (int i{1}; i < m_u.x_max() - 1; ++i)
			m_u(i, m_u.y_max() - 1) = 2 * bcTop[0] - m_u(i, m_u.y_max() - 2);
	else //< Neumann
		for (int i{1}; i < m_u.x_max() - 1; ++i)
		{
			m_u(i, m_u.y_max() - 1) = m_u(i, m_u.y_max() - 2);
			m_v(i, m_u.y_max() - 1) = m_v(i, m_u.y_max() - 2) + bcTop[1] * m_dy;
		}

	// left
	if (useDirichletBc[LEFT]) //< Dirichlet
		for (int j{0}; j < m_v.y_max(); ++j)
			m_v(0, j) = 2 * bcLeft[1] - m_v(1, j);
	else //< Neumann
		for (int j{0}; j < m_v.y_max(); ++j)
		{
			m_v(0, j) = m_v(1, j);
			m_u(0, j) = m_u(1, j) - bcLeft[0] * m_dx;
		}

	// right
	if (useDirichletBc[RIGHT]) //< Dirichlet
		for (int j{0}; j < m_v.y_max(); ++j)
			m_v(m_v.x_max() - 1, j) = 2 * bcRight[1] - m_v(m_v.x_max() - 2, j);
	else //< Neumann
		for (int j{0}; j < m_v.y_max(); ++j)
		{
			m_v(m_v.x_max() - 1, j) = m_v(m_v.x_max() - 2, j);
			m_u(m_v.x_max() - 1, j) = m_u(m_v.x_max() - 2, j) + bcRight[0] * m_dx;
		}

	/*
	// only the boundary values that depend on inner u,v are updated
	// update u
	for (int i{1}; i < m_u.x_max() - 1; ++i)
	{
		// bottom boundary
		m_u(i, 0) = 2 * bcBottom[0] - m_u(i, 1); //ok
		// top boundary
		m_u(i, m_u.y_max() - 1) = 2 * bcTop[0] - m_u(i, m_u.y_max() - 2);
	}

	// update v
	for (int j{0}; j < m_v.y_max(); ++j)
	{
		// left boundary
		m_v(0, j) = 2 * bcLeft[1] - m_v(1, j);
		// right boundary
		m_v(m_v.x_max() - 1, j) = 2 * bcRight[1] - m_v(m_v.x_max() - 2, j);
	}
	*/
}

void Discretization::compute_bound_val_FG()
{
	//! assumption: F^(n)_ij = u^(n+1)_ij and G^(n)_ij = v^(n+1)_ij
	//! set F
	for (int j{0}; j < m_u.y_max() - 1; ++j)
	{
		//! left boundary
		m_F(0, j) = m_u(0, j);
		//! right boundary
		m_F(m_u.x_max() - 1, j) = m_u(m_u.x_max() - 1, j);
	}
	//! set G
	for (int i{0}; i < m_v.x_max() - 1; ++i)
	{
		//! bottom boundary
		m_G(i, 0) = m_v(i, 0);
		//! top boundary
		m_G(i, m_v.y_max() - 1) = m_v(i, m_v.y_max() - 1);
	}
}

//! compute and update F and G
void Discretization::compute_FG()
{
	//! compute F
	for (int j{1}; j < m_u.y_max() - 1; ++j)
	{
		for (int i{1}; i < m_u.x_max() - 1; ++i)
		{
			m_F(i, j) = m_u(i, j)														   //< u at i,j
						+ m_dt * (1 / m_re * (compute_du_dx2(i, j) + compute_du_dy2(i, j)) //< diffusion term
								  - compute_du2_dx(i, j) - compute_duv_dy(i, j)			   //< convection terms
								  + m_g[0]);											   //< external force term (g_y)
		}
	}

	//! compute G
	for (int j{1}; j < m_v.y_max() - 1; ++j)
	{
		for (int i{1}; i < m_v.x_max() - 1; ++i)
		{
			m_G(i, j) = m_v(i, j)														   //< v at i,j
						+ m_dt * (1 / m_re * (compute_dv_dx2(i, j) + compute_dv_dy2(i, j)) //< diffusion term
								  - compute_duv_dx(i, j) - compute_dv2_dy(i, j)			   //< convection terms
								  + m_g[1]);											   //< external force term (g_x)
		}
	}
}

//! compute and update RHS
void Discretization::compute_RHS()
{
	assert(m_RHS.size() == m_F.size());
	assert(m_RHS.size() == m_G.size());

	for (int j{1}; j < m_RHS.size()[1] - 1; ++j)
	{
		for (int i{1}; i < m_RHS.size()[0] - 1; ++i)
		{
			m_RHS(i, j) = 1 / m_dt * ((m_F(i, j) - m_F(i - 1, j)) / m_dx	 //< momentum u
									  + (m_G(i, j) - m_G(i, j - 1)) / m_dy); //< momentum v
		}
	}
}

//! compute and update u and v
void Discretization::compute_uv()
{
	//! compute final u^(n+1)
	for (int j{1}; j < m_u.y_max() - 1; ++j)
	{
		for (int i{1}; i < m_u.x_max() - 1; ++i)
		{
			m_u(i, j) = m_F(i, j) - m_dt * compute_dp_dx(i, j);
		}
	}

	//! compute final v^(n+1)
	for (int j{1}; j < m_v.y_max() - 1; ++j)
	{
		for (int i{1}; i < m_v.x_max() - 1; ++i)
		{
			m_v(i, j) = m_G(i, j) - m_dt * compute_dp_dy(i, j);
		}
	}
}