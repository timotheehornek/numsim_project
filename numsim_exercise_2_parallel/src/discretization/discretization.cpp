
#include "discretization.h"

Discretization::Discretization(const std::array<int, 2> &nCells_global, const std::array<double, 2> &physicalSize, const double re, const std::array<double, 2> &g, Message_passer_2D &MP)
	: m_ex_nbrs{
		  ((MP.size() - MP.prcs(0)) <= MP.rank() ? false : true),		 //< top
		  ((MP.rank() % MP.prcs(0)) == 0 ? false : true),				 //< left
		  ((MP.rank() % MP.prcs(0)) == (MP.prcs(0) - 1) ? false : true), //< right
		  (MP.rank() < MP.prcs(0) ? false : true)						 //< bottom
	  },

	  m_nbrs{
		  (ex_nbrs(SIDE_TOP) ? (MP.rank() + MP.prcs(0)) : -1),
		  (ex_nbrs(SIDE_LEFT) ? (MP.rank() - 1) : -1),
		  (ex_nbrs(SIDE_RIGHT) ? (MP.rank() + 1) : -1),
		  (ex_nbrs(SIDE_BOTTOM) ? (MP.rank() - MP.prcs(0)) : -1),
	  },

	  m_nCells{
		  ex_nbrs(SIDE_RIGHT) ? nCells_global[0] / MP.prcs(0) : nCells_global[0] / MP.prcs(0) + nCells_global[0] % MP.prcs(0),
		  ex_nbrs(SIDE_TOP) ? nCells_global[1] / MP.prcs(1) : nCells_global[1] / MP.prcs(1) + nCells_global[1] % MP.prcs(1),
	  },

	  m_u{m_nCells[0] + 2, m_nCells[1] + 2}, m_v{m_nCells[0] + 2, m_nCells[1] + 2}, m_p{m_nCells[0] + 2, m_nCells[1] + 2}, m_F{m_nCells[0] + 2, m_nCells[1] + 2}, m_G{m_nCells[0] + 2, m_nCells[1] + 2}, m_RHS{m_nCells[0] + 2, m_nCells[1] + 2},

	  m_dx{physicalSize[0] / nCells_global[0]}, m_dy{physicalSize[1] / nCells_global[1]},

	  m_re{re},

	  m_g{g},

	  m_MP{MP}
{
}

//! compute dt from stability conditions
void Discretization::set_dt(double dt_max, double sec_factor)
{
	//! gather all abs max of u and v in all processes
	const double u_abs_max_local{m_u.abs_max()};
	const double v_abs_max_local{m_v.abs_max()};

	double u_abs_max_global;
	double v_abs_max_global;
	m_MP.reduce_max(&u_abs_max_global, u_abs_max_local);
	m_MP.reduce_max(&v_abs_max_global, v_abs_max_local);

	//! compute dt in root
	if (m_MP.rank() == 0)
	{
		//! stability condition induced by diffusion operator (identical in all processes)
		const double stab_cond_1{
			m_re / 2 *
			(std::pow(m_dx, 2) * std::pow(m_dy, 2)) / (std::pow(m_dx, 2) + std::pow(m_dy, 2))};

		//! stability conditions induced by convection operator
		const double stab_cond_2{m_dx / u_abs_max_global};
		const double stab_cond_3{m_dy / v_abs_max_global};

		//! find largest possible dt considering security factor
		const double dt_best{std::min(std::min(stab_cond_1, stab_cond_2), stab_cond_3) * sec_factor};

		//! set dt respecting user input for highest dt
		m_dt = std::min(dt_best, dt_max);
	}
	//! broadcast dt from process 0 to other processes
	m_MP.broadcast(&m_dt);

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
bool Discretization::ex_nbrs(int dir) const
{
	assert(0 <= dir && dir <= 3);
	return m_ex_nbrs[dir];
}

//! returns ranks of neighbors (-1 if not existent)
//! (directions: 0=top, 1=left, 2=right, 3=bottom)
int Discretization::nbrs(int dir) const
{
	return m_nbrs[dir];
}

void Discretization::communicate_update_bound_val_uv(
	const std::array<double, 2> &dirichletBcBottom, const std::array<double, 2> &dirichletBcTop,
	const std::array<double, 2> &dirichletBcLeft, const std::array<double, 2> &dirichletBcRight)
{
	//! declare receive buffers
	std::vector<double> u_top_buffer;
	std::vector<double> u_bottom_buffer;
	std::vector<double> u_left_buffer;
	std::vector<double> u_right_buffer;

	std::vector<double> v_top_buffer;
	std::vector<double> v_bottom_buffer;
	std::vector<double> v_left_buffer;
	std::vector<double> v_right_buffer;

	double u_top_left_buffer;
	double v_bottom_right_buffer;

	//! send boundaries or set boundaries
	//! top
	if (ex_nbrs(SIDE_TOP))
	{
		//! left neighbour exists -> send u and v to neighbour

		//! extract and send u
		m_MP.send_top(m_u.get_side(SIDE_TOP), TAG_U);

		//! receive u
		u_top_buffer.resize(m_nCells[0]);
		m_MP.receive_top(u_top_buffer, TAG_U);

		//! extract and send v
		m_MP.send_top(m_v.get_side(SIDE_TOP), TAG_V);

		//! receive v
		v_top_buffer.resize(m_nCells[0]);
		m_MP.receive_top(v_top_buffer, TAG_V);
	}

	//! bottom
	if (ex_nbrs(SIDE_BOTTOM))
	{
		//! left neighbour exists -> send u and v to neighbour

		//! extract and send u
		m_MP.send_bottom(m_u.get_side(SIDE_BOTTOM), TAG_U);

		//! receive u
		u_bottom_buffer.resize(m_nCells[0]);
		m_MP.receive_bottom(u_bottom_buffer, TAG_U);

		//! extract and send v
		m_MP.send_bottom(m_v.get_side(SIDE_BOTTOM), TAG_V);

		//! receive v
		v_bottom_buffer.resize(m_nCells[0]);
		m_MP.receive_bottom(v_bottom_buffer, TAG_V);
	}

	//! right
	if (ex_nbrs(SIDE_RIGHT))
	{
		//! right neighbour exists -> send u and v to neighbour

		//! extract and send u
		m_MP.send_right(m_u.get_side(SIDE_RIGHT), TAG_U);

		//! receive u
		u_right_buffer.resize(m_nCells[1]);
		m_MP.receive_right(u_right_buffer, TAG_U);

		//! extract and send v
		m_MP.send_right(m_v.get_side(SIDE_RIGHT), TAG_V);

		//! receive v
		v_right_buffer.resize(m_nCells[1]);
		m_MP.receive_right(v_right_buffer, TAG_V);
	}

	//! left
	if (ex_nbrs(SIDE_LEFT))
	{
		//! left neighbour exists -> send u and v to neighbour

		//! extract and send u
		m_MP.send_left(m_u.get_side(SIDE_LEFT), TAG_U);

		//! receive u
		u_left_buffer.resize(m_nCells[1]);
		m_MP.receive_left(u_left_buffer, TAG_U);

		//! extract and send v
		m_MP.send_left(m_v.get_side(SIDE_LEFT), TAG_V);

		//! receive v
		v_left_buffer.resize(m_nCells[1]);
		m_MP.receive_left(v_left_buffer, TAG_V);
	}

	//! corners
	if (ex_nbrs(SIDE_LEFT) && ex_nbrs(SIDE_TOP))
	{
		m_MP.send_top_left(m_v.get_corner(CORNER_TOP_LEFT), TAG_V);

		m_MP.receive_top_left(u_top_left_buffer, TAG_U);
	}
	if (ex_nbrs(SIDE_RIGHT) && ex_nbrs(SIDE_BOTTOM))
	{
		m_MP.send_bottom_right(m_u.get_corner(CORNER_BOTTOM_RIGHT), TAG_U);

		m_MP.receive_bottom_right(v_bottom_right_buffer, TAG_V);
	}

	//! wait for message passing to complete
	//m_MP.wait();

	//! write boundaries
	//! right
	if (ex_nbrs(SIDE_RIGHT))
	{
		m_u.write_bound(SIDE_RIGHT, u_right_buffer);
		m_v.write_bound(SIDE_RIGHT, v_right_buffer);
	}

	//! left
	if (ex_nbrs(SIDE_LEFT))
	{
		m_u.write_bound(SIDE_LEFT, u_left_buffer);
		m_v.write_bound(SIDE_LEFT, v_left_buffer);
	}

	//! top
	if (ex_nbrs(SIDE_TOP))
	{
		m_u.write_bound(SIDE_TOP, u_top_buffer);
		m_v.write_bound(SIDE_TOP, v_top_buffer);
	}

	//! bottom
	if (ex_nbrs(SIDE_BOTTOM))
	{
		m_u.write_bound(SIDE_BOTTOM, u_bottom_buffer);
		m_v.write_bound(SIDE_BOTTOM, v_bottom_buffer);
	}

	//! corners
	if (ex_nbrs(SIDE_LEFT) && ex_nbrs(SIDE_TOP))
	{
		m_u.write_corner(CORNER_TOP_LEFT, u_top_left_buffer);
	}
	if (ex_nbrs(SIDE_RIGHT) && ex_nbrs(SIDE_BOTTOM))
	{
		m_v.write_corner(CORNER_BOTTOM_RIGHT, v_bottom_right_buffer);
	}

	//! compute boundaries using boundary conditions
	if (!ex_nbrs(SIDE_TOP))
	{
		//! do normal boundary treatment for u,v
		//! top boundary
		//! set u
		for (int i{0}; i < m_u.x_max(); ++i)
		{
			m_u(i, m_u.y_max() - 1) = 2 * dirichletBcTop[0] - m_u(i, m_u.y_max() - 2);
		}
		//! set v
		for (int i{0}; i < m_v.x_max(); ++i)
		{
			m_v(i, m_v.y_max() - 2) = dirichletBcTop[1];
		}
	}

	//! bottom
	if (!ex_nbrs(SIDE_BOTTOM))
	{
		//! do normal boundary treatment for u,v
		//! bottom boundary
		//! set u
		for (int i{0}; i < m_u.x_max(); ++i)
		{
			m_u(i, 0) = 2 * dirichletBcBottom[0] - m_u(i, 1);
		}
		//! set v
		for (int i{0}; i < m_v.x_max(); ++i)
		{
			m_v(i, 0) = dirichletBcBottom[1];
		}
	}
	//! right
	if (!ex_nbrs(SIDE_RIGHT))
	{
		//! do normal boundary treatment for u,v
		//! right boundary
		//! set u
		for (int j{0}; j < m_u.y_max(); ++j) //end one above
		{
			m_u(m_u.x_max() - 2, j) = dirichletBcRight[0];
		}

		//! set v
		for (int j{0}; j < m_v.y_max(); ++j) //start one below
		{
			m_v(m_v.x_max() - 1, j) = 2 * dirichletBcRight[1] - m_v(m_v.x_max() - 2, j);
		}
	}

	//! left
	if (!ex_nbrs(SIDE_LEFT))
	{
		//! do normal boundary treatment for u,v
		//! left boundary
		//! set u
		for (int j{0}; j < m_u.y_max(); ++j) //start one below and end one above
		{
			m_u(0, j) = dirichletBcLeft[0];
		}
		//! set v
		for (int j{0}; j < m_v.y_max(); ++j) //start one below
		{
			m_v(0, j) = 2 * dirichletBcLeft[1] - m_v(1, j);
		}
	}
}

//! compute and update F and G
void Discretization::compute_FG()
{
	//! compute F
	if (!ex_nbrs(SIDE_RIGHT))
	{
		//! do not compute the most right coloumn
		for (int j{1}; j < m_u.y_max() - 1; ++j)
		{
			for (int i{1}; i < m_u.x_max() - 2; ++i) //< one less
			{
				m_F(i, j) = m_u(i, j)														   //< u at i,j
							+ m_dt * (1 / m_re * (compute_du_dx2(i, j) + compute_du_dy2(i, j)) //< diffusion term
									  - compute_du2_dx(i, j) - compute_duv_dy(i, j)			   //< convection terms
									  + m_g[0]);											   //< external force term (g_y)
			}
		}
	}
	else
	{
		//! compute F for all inner cells
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
	}

	//! compute G
	if (!ex_nbrs(SIDE_TOP))
	{
		//! do not compute the most top row
		for (int j{1}; j < m_v.y_max() - 2; ++j) //< one less
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
	else
	{
		//! compute G for all inner cells
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
}

//! communicate_FG
void Discretization::communicate_update_bound_val_FG()
{
	//! declare receive buffers
	std::vector<double> F_left_buffer;
	std::vector<double> G_bottom_buffer;

	//! send F if right neighbour exists
	if (ex_nbrs(SIDE_RIGHT))
	{
		m_MP.send_right(m_F.get_side(SIDE_RIGHT), TAG_F);
	}
	else
	{
		for (int j{1}; j < m_u.y_max() - 1; ++j)
		{
			//! right boundary
			m_F(m_u.x_max() - 2, j) = m_u(m_u.x_max() - 2, j);
		}
	}

	//! send G if top neighbour exists
	if (ex_nbrs(SIDE_TOP))
	{
		m_MP.send_top(m_G.get_side(SIDE_TOP), TAG_G);
	}
	else
	{
		for (int i{1}; i < m_v.x_max() - 1; ++i)
		{
			//! top boundary
			m_G(i, m_v.y_max() - 2) = m_v(i, m_v.y_max() - 2);
		}
	}

	//! receive if left neighbour exists
	if (ex_nbrs(SIDE_LEFT))
	{
		F_left_buffer.resize(m_nCells[1]);
		m_MP.receive_left(F_left_buffer, TAG_F);
	}
	else
	{
		for (int j{1}; j < m_u.y_max() - 1; ++j)
		{
			//! left boundary
			m_F(0, j) = m_u(0, j);
		}
	}
	//! receive if bottom neighbour exists
	if (ex_nbrs(SIDE_BOTTOM))
	{
		G_bottom_buffer.resize(m_nCells[0]);
		m_MP.receive_bottom(G_bottom_buffer, TAG_G);
	}
	else
	{
		for (int i{1}; i < m_v.x_max() - 1; ++i)
		{
			//! bottom boundary
			m_G(i, 0) = m_v(i, 0);
		}
	}

	//! wait for message passing to complete
	//m_MP.wait();

	//! write boundaries
	if (ex_nbrs(SIDE_LEFT))
	{
		m_F.write_bound(SIDE_LEFT, F_left_buffer);
	}
	if (ex_nbrs(SIDE_BOTTOM))
	{
		m_G.write_bound(SIDE_BOTTOM, G_bottom_buffer);
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