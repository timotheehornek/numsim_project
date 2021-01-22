#include "lattice_boltzmann.h"
/*
Directions:
0 - center
1 - right
2 - top
3 - left
4 - bottom
5 - top-right
6 - top-left
7 - bottom-left
8 - bottom-right
*/

Lattice_boltzmann::Lattice_boltzmann(const std::array<int, 2>& nCells, const std::array<double, 2>& physicalSize, const double re, const std::array<double, 2>& g)
	: m_nCells{ nCells },
	m_u{ nCells[0], nCells[1]},
	m_v{ nCells[0], nCells[1]},
	m_p{ nCells[0], nCells[1]},

	m_obstacle{ nCells[0], nCells[1]},

	m_dx{ physicalSize[0] / (nCells[0]) },
	m_dy{ physicalSize[1] / (nCells[1]) },
	m_re{ re },
	m_g{ g },

	m_lattice_factor { 1.0 }
{}

//! set dt by taking dt from user
void Lattice_boltzmann::set_dt(double dt_max)
{
	//! set dt respecting user input for highest dt
	m_dt = dt_max;

	assert(m_dt > 0);
}

//! set dt according to lattice parameter
void Lattice_boltzmann::set_dt()
{
	assert(m_dx == m_dy);

	//set dt according to lattice factor
	m_dt = m_dx * m_lattice_factor;

	assert(m_dt > 0);
}

//! compute the relaxtion parameter tau
void Lattice_boltzmann::set_tau(double viscosity)
{
	//compute tau
	double tau = 0.5 + 3.0 * viscosity / m_dx;
	//set one over tau
	m_one_over_tau = 1.0 / tau;

	//print relaxation factor
	std::cout << "1/tau = " << m_one_over_tau << '\n';
	assert(m_one_over_tau < 2);
}

//! manually set the relaxtion parameter tau
void Lattice_boltzmann::set_tau_manually(double tau)
{
	//set one over tau
	m_one_over_tau = 1.0/ tau;

	assert(m_one_over_tau < 2);
}

//! compute the viscosity of the fluid - relevant for setting tau
double Lattice_boltzmann::compute_viscosity(double length, const std::array<bool, 4>& useDirichletBc,
		const std::array<double, 2>& bcBottom, const std::array<double, 2>& bcTop,
		const std::array<double, 2>& bcLeft, const std::array<double, 2>& bcRight)
{
	double viscosity = 0;

	//first check wheater pressure inflow is used
	// -> use pressure gradient to determine viscosity
	if(!useDirichletBc[0] && !useDirichletBc[1] && std::abs(bcBottom[1] - bcTop[1]) > 0) // bottom(0) and top(1)
	{
		viscosity = std::sqrt(length * length * length * std::abs(bcBottom[1] - bcTop[1]) / (8 * m_re));
	}
	else if (!useDirichletBc[2] && !useDirichletBc[3] && std::abs(bcLeft[0] - bcRight[0]) > 0) // left(2) and right(3)
	{
		viscosity = std::sqrt(length * length * length * std::abs(bcLeft[0] - bcRight[0]) / (8 * m_re));
	}

	//as there are only velocity inflow or no-slip boundaries use
	// the maximal velocity to compute the viscosity
	else
	{
		//compute maximal velocity on boundary
		double max_vel = get_max_vel(useDirichletBc, bcBottom, bcTop, bcLeft, bcRight);
		//check if maximal velocity is creater than 0
		if (max_vel > 0)
		{
			viscosity =  max_vel * length / m_re;
		}
		// else set viscosity such that relaxation parameter tau = 1;
		else
		{
			std::cout << "viscosity could no be computed -> tau set to 1" << '\n';
			viscosity = m_dx / 6.0;
		}
	}

	//print viscosity
	std::cout << "viscosity = " << viscosity << '\n';
	return viscosity;
}

//! set m_obstacle to match a rectangular object
void Lattice_boltzmann::set_obstacle_rect(const std::array<int, 4>& obstaclePos)
{
	for (int i = obstaclePos[0]; i < obstaclePos[2]; i++)
	{
		for (int j = obstaclePos[1]; j < obstaclePos[3]; j++)
		{
			// set obstacle bool on true
			m_obstacle(i,j) = 1.0;
		}
	}
}

//! set m_obstacle to match a rectangular object
void Lattice_boltzmann::set_obstacle_round(const std::array<int, 2>& obstacleCenter, double radius)
{
	for (int i = obstacleCenter[0] - radius ; i <= obstacleCenter[0] + radius; i++)
	{
		for (int j = obstacleCenter[1] - radius; j <= obstacleCenter[1] + radius; j++)
		{
			double distance = std::sqrt((obstacleCenter[0] - i) * (obstacleCenter[0] - i)
											+ (obstacleCenter[1] - j) * (obstacleCenter[1] - j));
			// if in radius set obstacle bool on true
			if (distance < radius)
			{
				m_obstacle(i,j) = 1.0;
			}
		}
	}
}

//! get max velocity set on boundary
double Lattice_boltzmann::get_max_vel(const std::array<bool, 4>& useDirichletBc,
		const std::array<double, 2>& bcBottom, const std::array<double, 2>& bcTop,
		const std::array<double, 2>& bcLeft, const std::array<double, 2>& bcRight) const
{
	double max_vel = 0.0;
	//check bottom boundary
	if (useDirichletBc[0])
	{
		if (std::abs(bcBottom[0]) > max_vel)
		{
			max_vel = std::abs(bcBottom[0]);
		}
		if (std::abs(bcBottom[1]) > max_vel)
		{
			max_vel = std::abs(bcBottom[1]);
		}
	}
	//check top boundary
	if (useDirichletBc[1])
	{
		if (std::abs(bcTop[0]) > max_vel)
		{
			max_vel = std::abs(bcTop[0]);
		}
		if (std::abs(bcTop[1]) > max_vel)
		{
			max_vel = std::abs(bcTop[1]);
		}
	}
	//check left boundary
	if (useDirichletBc[2])
	{
		if (std::abs(bcLeft[0]) > max_vel)
		{
			max_vel = std::abs(bcLeft[0]);
		}
		if (std::abs(bcLeft[1]) > max_vel)
		{
			max_vel = std::abs(bcLeft[1]);
		}
	}
	//check right boundary
	if (useDirichletBc[3])
	{
		if (std::abs(bcRight[0]) > max_vel)
		{
			max_vel = std::abs(bcRight[0]);
		}
		if (std::abs(bcRight[1]) > max_vel)
		{
			max_vel = std::abs(bcRight[1]);
		}
	}

	return max_vel;
}

//! initialize f_vector in equilibrium state
void Lattice_boltzmann::initialize_f(double density, std::array<Array2D, 9>& input)
{
	for(int i = 0; i < m_nCells[0]; i++)
	{
		for(int j = 0; j < m_nCells[1]; j++)
		{
			input[0](i,j) = 4.0 / 9.0 * density;

			input[1](i,j) = 1.0 / 9.0 * density;
			input[2](i,j) = 1.0 / 9.0 * density;
			input[3](i,j) = 1.0 / 9.0 * density;
			input[4](i,j) = 1.0 / 9.0 * density;

			input[5](i,j) = 1.0 / 36.0 * density;
			input[6](i,j) = 1.0 / 36.0 * density;
			input[7](i,j) = 1.0 / 36.0 * density;
			input[8](i,j) = 1.0 / 36.0 * density;
		}
	}
}

//! getter
double Lattice_boltzmann::dt() const
{
	return m_dt;
}
double Lattice_boltzmann::dx() const
{
	return m_dx;
}
double Lattice_boltzmann::dy() const
{
	return m_dy;
}
const Staggered_grid& Lattice_boltzmann::u() const
{
	return m_u;
}
const Staggered_grid& Lattice_boltzmann::v() const
{
	return m_v;
}
const Staggered_grid& Lattice_boltzmann::p() const
{
	return m_p;
}
const std::array<int, 2>& Lattice_boltzmann::nCells() const
{
	return m_nCells;
}
const double Lattice_boltzmann::u(int i, int j) const
{
	return m_u(i,j);
}
const double Lattice_boltzmann::v(int i, int j) const
{
	return m_v(i,j);
}
const double Lattice_boltzmann::p(int i, int j) const
{
	return m_p(i,j);
}

//! compute rho by summing over all f
double Lattice_boltzmann::compute_rho(int i, int j, const std::array<Array2D, 9>& input) const
{
  double rho = 0;
  for (int k = 0; k < 9; k++)
  {
    rho += input[k](i,j);
  }
  return rho;
}

//! compute u by multipliing all f with the first components of the corresponding lattice vectors
double Lattice_boltzmann::compute_u(int i, int j, const std::array<Array2D, 9>& input) const
{
  double vel_u = input[1](i,j) + input[5](i,j) + input[8](i,j)
               - (input[3](i,j) + input[6](i,j) + input[7](i,j));
  return vel_u;
}

//! compute v by multipliing all f with the first components of the corresponding lattice vectors
double Lattice_boltzmann::compute_v(int i, int j, const std::array<Array2D, 9>& input) const
{
  double vel_v = input[2](i,j) + input[5](i,j) + input[6](i,j)
           		 - (input[4](i,j) + input[7](i,j) + input[8](i,j));
  return vel_v;
}

//! compute the f_eq and the relaxate towards the equilibrium
void Lattice_boltzmann::collide(int i, int j, double density, double vel_u, double vel_v, std::array<Array2D, 9>& input)
{
	//stores values for visualization
	m_p(i,j) = density; //density / (3 * m_lattice_factor * m_lattice_factor);
	m_u(i,j) = vel_u;
	m_v(i,j) = vel_v;

	//compute equilibrium values
	double f_eq;

  //compute absolute velocity term
  double abs_vel = 1.5f * m_lattice_factor * m_lattice_factor * (vel_u * vel_u + vel_v * vel_v);

  //center
  //compute equilibrium
  f_eq = 4.0 / 9.0 * (density - abs_vel);
  //relaxte towards equilibrium
	input[0](i,j) = (1.0 - m_one_over_tau) * input[0](i,j) + m_one_over_tau * f_eq;

  //right
  //compute equilibrium
  f_eq = 1.0 / 9.0 * (density + 3.0 * m_lattice_factor * vel_u + 4.5 * m_lattice_factor * m_lattice_factor * vel_u * vel_u - abs_vel);
  //relaxte towards equilibrium
	input[1](i,j) = (1.0 - m_one_over_tau) * input[1](i,j) + m_one_over_tau * f_eq;

  //top
  //compute equilibrium
  f_eq = 1.0 / 9.0 * (density + 3.0 * m_lattice_factor * vel_v + 4.5 * m_lattice_factor * m_lattice_factor * vel_v * vel_v - abs_vel);
  //relaxte towards equilibrium
	input[2](i,j) = (1.0 - m_one_over_tau) * input[2](i,j) + m_one_over_tau * f_eq;

  //left
  //compute equilibrium
  f_eq = 1.0 / 9.0 * (density - 3.0 * m_lattice_factor * vel_u + 4.5 * m_lattice_factor * m_lattice_factor * vel_u * vel_u - abs_vel);
  //relaxte towards equilibrium
	input[3](i,j) = (1.0 - m_one_over_tau) * input[3](i,j) + m_one_over_tau * f_eq;

  //bottom
  //compute equilibrium
  f_eq = 1.0 / 9.0 * (density - 3.0 * m_lattice_factor * vel_v + 4.5 * m_lattice_factor * m_lattice_factor * vel_v * vel_v - abs_vel);
  //relaxte towards equilibrium
	input[4](i,j) = (1.0 - m_one_over_tau) * input[4](i,j) + m_one_over_tau * f_eq;

  //top-right
  //compute equilibrium
  f_eq = 1.0 / 36.0 * (density + 3.0 * m_lattice_factor * (vel_u + vel_v) + 4.5 * m_lattice_factor * m_lattice_factor * (vel_u + vel_v) * (vel_u + vel_v) - abs_vel);
  //relaxte towards equilibrium
	input[5](i,j) = (1.0 - m_one_over_tau) * input[5](i,j) + m_one_over_tau * f_eq;

  //top-left
  //compute equilibrium
  f_eq = 1.0 / 36.0 * (density + 3.0 * m_lattice_factor * (-vel_u + vel_v) + 4.5 * m_lattice_factor * m_lattice_factor * (-vel_u + vel_v) * (-vel_u + vel_v) - abs_vel);
	//relaxte towards equilibrium
	input[6](i,j) = (1.0 - m_one_over_tau) * input[6](i,j) + m_one_over_tau * f_eq;

  //bottom-left
  //compute equilibrium
  f_eq = 1.0 / 36.0 * (density - 3.0 * m_lattice_factor * (vel_u + vel_v) + 4.5 * m_lattice_factor * m_lattice_factor * (vel_u + vel_v) * (vel_u + vel_v) - abs_vel);
	//relaxte towards equilibrium
	input[7](i,j) = (1.0 - m_one_over_tau) * input[7](i,j) + m_one_over_tau * f_eq;

  //bottom right
  //compute equilibrium
  f_eq = 1.0 / 36.0 * (density + 3.0 * m_lattice_factor * (vel_u - vel_v) + 4.5 * m_lattice_factor * m_lattice_factor * (vel_u - vel_v) * (vel_u - vel_v) - abs_vel);
  //relaxte towards equilibrium
	input[8](i,j) = (1.0 - m_one_over_tau) * input[8](i,j) + m_one_over_tau * f_eq;
}

//! bounceback (swap) values as object cell at i,j
void Lattice_boltzmann::bounceback(int i, int j, std::array<Array2D, 9>& input)
{
  //full-way bounceback opposite values
  double holder;

  //swap left and right
  holder = input[1](i,j);
  input[1](i,j) = input[3](i,j);
  input[3](i,j) = holder;

  //swap top and bottom
  holder = input[2](i,j);
  input[2](i,j) = input[4](i,j);
  input[4](i,j) = holder;

  //swap top-right and bottom-left
  holder = input[5](i,j);
  input[5](i,j) = input[7](i,j);
  input[7](i,j) = holder;

  //swap top-left and bottom-right
  holder = input[6](i,j);
  input[6](i,j) = input[8](i,j);
  input[8](i,j) = holder;
}

void Lattice_boltzmann::stream(int i, int j, const std::array<Array2D, 9>& input, std::array<Array2D, 9>& output) const
{
  //send input values to neighbor cells in output vector
  output[0](i, j) = input[0](i, j);

  //send right value to right neighbor
  output[1](i + 1, j) = input[1](i, j);

  //send top value to top neighbor
  output[2](i, j + 1) = input[2](i, j);

	//send left value to left neighbor
	output[3](i - 1, j) = input[3](i, j);

	//send bottom value to bottom neighbor
	output[4](i, j - 1) = input[4](i, j);

	//send top-right value to top-right neighbor
	output[5](i + 1, j + 1) = input[5](i, j);

	//send top-left value to top-left neighbor
	output[6](i - 1, j + 1) = input[6](i, j);

	//send bottom-left value to bottom-left neighbor
	output[7](i - 1, j - 1) = input[7](i, j);

	//send bottom-right value to bottom-right neighbor
	output[8](i + 1, j - 1) = input[8](i, j);
}

//! run two lattice iteration with input and output on all inner nodes
void Lattice_boltzmann::two_lattice(std::array<Array2D, 9>& input, std::array<Array2D, 9>& output)
{
  //get loop size
  int size_x = input[0].size()[0];
  int size_y = input[0].size()[1];

  //loop over all inner nodes
  for(int i = 1; i < size_x - 1; i++)
  {
    for(int j = 1; j < size_y - 1; j++)
    {
			//check if obstacle and choose bounceback or collide step
      if(m_obstacle(i,j))
      {
        //obstacle cells bounceback instead collide step
        bounceback(i, j, input);
      }
      else
      {
        //collide step
        collide(i, j, compute_rho(i,j, input), compute_u(i,j, input), compute_v(i,j, input), input);
      }

			//stream step
			stream(i, j, input, output);
		}
	}
}

void Lattice_boltzmann::boundary_treatment(
	std::array<Array2D, 9>& input, std::array<Array2D, 9>& output,
	const std::array<bool, 4>& useDirichletBc,
	const std::array<double, 2>& bcBottom, const std::array<double, 2>& bcTop,
	const std::array<double, 2>& bcLeft, const std::array<double, 2>& bcRight)
{
  //General boundary conditions

	//std::array<bool, 4> useDirichletBc; //< {bottom, top, left, right}
	//< type of boundary condition true for Dirichlet and false for Neumann

  //get loop size
  int size_x = input[0].size()[0];
  int size_y = input[0].size()[1];

	/////////////////////////////////////////////////////////////////////////////
	//CORNERS
	double i, j, density, f_buried;
  //////////////////////////////////////////////////////////////////////////////
  //top-right corner
	i = size_x - 1;
	j = size_y - 1;
	//compute missing values
	//half-way bounceback the values of the inner cell pointing out of the domain
	input[3](i, j) = input[1](i, j);
	input[4](i, j) = input[2](i, j);
	input[7](i, j) = input[5](i, j);

	//compute buried link
	f_buried = 1.0 / 18.0 * (input[1](i, j) + input[2](i, j) + input[3](i, j) + input[4](i, j) + input[5](i, j) + input[7](i, j));

	if (!useDirichletBc[3])
	{//outflow boundary set values according to set density
		//initialize as outflow
		double density = 1.0;
		//check if inflow - assume full pressure gradient is applied to inflow
		if (bcRight[0] > bcLeft[0])
		{
			density += 3 * (bcRight[0] - bcLeft[0]) / m_lattice_factor * size_x;
		}

		f_buried = density / 18.0 - f_buried;
	}
	else if (!useDirichletBc[1])
	{//outflow boundary set values according to set density
		//initialize as outflow
		double density = 1.0;
		//check if inflow - assume full pressure gradient is applied to inflow
		if (bcTop[1] > bcBottom[1])
		{
			density += 3 * (bcTop[1] - bcBottom[1]) / m_lattice_factor * size_y;
		}

		f_buried = density / 18.0 - f_buried;
	}

	input[6](i, j) = f_buried;
	input[8](i, j) = f_buried;
	input[0](i, j)= 16.0 * f_buried;

	//normal collide step
	collide(i, j, compute_rho(i,j, input), compute_u(i,j, input), compute_v(i,j, input), input);

	//modified stream step -> send values into domain
	//send input values to neighbor cells in output vector
	output[3](i - 1, j) = input[3](i, j);
	output[4](i, j - 1) = input[4](i, j);
	output[7](i - 1, j - 1) = input[7](i, j);

	//////////////////////////////////////////////////////////////////////////////
  //top-left corner
	i = 0;
	j = size_y - 1;

	//compute missing values
	//half-way bounceback the values of the inner cell pointing out of the domain
	input[1](i, j) = input[3](i, j);
	input[4](i, j) = input[2](i, j);
	input[8](i, j) = input[6](i, j);

	//compute buried link
	f_buried = 1.0 / 18.0 * (input[1](i, j) + input[2](i, j) + input[3](i, j) + input[4](i, j) + input[6](i, j) + input[8](i, j));

	if (!useDirichletBc[2])
	{//outflow boundary set values according to set density
		//initialize as outflow
		double density = 1.0;
		//check if inflow - assume full pressure gradient is applied to inflow
		if (bcLeft[0] > bcRight[0])
		{
			density += 3 * (bcLeft[0] - bcRight[0]) / m_lattice_factor * size_x;
		}

	  f_buried = density / 18.0 - f_buried;
	}
	else if (!useDirichletBc[1])
	{//outflow boundary set values according to set density
		//initialize as outflow
		double density = 1.0;
		//check if inflow - assume full pressure gradient is applied to inflow
		if (bcTop[1] > bcBottom[1])
		{
			density += 3 * (bcTop[1] - bcBottom[1]) / m_lattice_factor * size_y;
		}

		f_buried = density / 18.0 - f_buried;
	}

	input[5](i, j) = f_buried;
	input[7](i, j) = f_buried;
	input[0](i, j)= 16.0 * f_buried;

	//normal collide step
	collide(i, j, compute_rho(i,j, input), compute_u(i,j, input), compute_v(i,j, input), input);

	//modified stream step -> send values into domain
	//send input values to neighbor cells in output vector
	output[1](i + 1, j) = input[1](i, j);
	output[4](i, j - 1) = input[4](i, j);
	output[8](i + 1, j - 1) = input[8](i, j);

	//////////////////////////////////////////////////////////////////////////////
  //bottom-left corner -> index (0,0)
	i = 0;
	j = 0;
	//compute missing values
	//half-way bounceback the values of the inner cell pointing out of the domain
	input[1](i, j) = input[3](i, j);
	input[2](i, j) = input[4](i, j);
	input[5](i, j) = input[7](i, j);

	//compute buried link
	f_buried = 1.0 / 18.0 * (input[1](i, j) + input[2](i, j) + input[3](i, j) + input[4](i, j) + input[5](i, j) + input[7](i, j));

	if (!useDirichletBc[2])
	{//outflow boundary set values according to set density
		//initialize as outflow
		double density = 1.0;
		//check if inflow - assume full pressure gradient is applied to inflow
		if (bcLeft[0] > bcRight[0])
		{
			density += 3 * (bcLeft[0] - bcRight[0]) / m_lattice_factor * size_x;
		}

	  f_buried = density / 18.0 - f_buried;
	}
	else if (!useDirichletBc[0])
	{//outflow boundary set values according to set density
		//initialize as outflow
		double density = 1.0;
		//check if inflow - assume full pressure gradient is applied to inflow
		if (bcBottom[1] > bcTop[1])
		{
			density += 3 * (bcBottom[1] - bcTop[1]) / m_lattice_factor * size_y;
		}


		f_buried = density / 18.0 - f_buried;
	}

	input[6](i, j) = f_buried;
	input[8](i, j) = f_buried;
	input[0](i,j)= 16.0 * f_buried;

	//normal collide step
	collide(i, j, compute_rho(i,j, input), compute_u(i,j, input), compute_v(i,j, input), input);

	//modified stream step -> send values into domain
	//send input values to neighbor cells in output vector
	output[1](i + 1, j) = input[1](i, j);
	output[2](i, j + 1) = input[2](i, j);
	output[5](i + 1, j + 1) = input[5](i, j);

	//////////////////////////////////////////////////////////////////////////////
  //bottom-right corner -> index (size_x - 1, 0)
	i = size_x - 1;
	j = 0;

	//compute missing values
	//half-way bounceback the values of the inner cell pointing out of the domain
	input[3](i, j) = input[1](i, j);
	input[2](i, j) = input[4](i, j);
	input[6](i, j) = input[8](i, j);

	//compute buried link
	f_buried = 1.0 / 18.0 * (input[1](i, j) + input[2](i, j) + input[3](i, j) + input[4](i, j) + input[6](i, j) + input[8](i, j));

	if (!useDirichletBc[3])
	{//outflow boundary set values according to set density
		//initialize as outflow
		double density = 1.0;
		//check if inflow - assume full pressure gradient is applied to inflow
		if (bcRight[0] > bcLeft[0])
		{
			density += 3 * (bcRight[0] - bcLeft[0]) / m_lattice_factor * size_x;
		}

		f_buried = density / 18.0 - f_buried;
	}
	else if (!useDirichletBc[0])
	{//outflow boundary set values according to set density
		//initialize as outflow
		double density = 1.0;
		//check if inflow - assume full pressure gradient is applied to inflow
		if (bcBottom[1] > bcTop[1])
		{
			density += 3 * (bcBottom[1] - bcTop[1]) / m_lattice_factor * size_y;
		}

		f_buried = density / 18.0 - f_buried;
	}

	input[5](i, j) = f_buried;
	input[7](i, j) = f_buried;
	input[0](i,j)= 16.0 * f_buried;

	//normal collide step
	collide(i, j, compute_rho(i,j, input), compute_u(i,j, input), compute_v(i,j, input), input);

	//modified stream step -> send values into domain
	//send input values to neighbor cells in output vector
	output[2](i, j + 1) = input[2](i, j);
	output[3](i - 1, j) = input[3](i, j);
	output[6](i - 1, j + 1) = input[6](i, j);

	/////////////////////////////////////////////////////////////////////////////
	//EDGES
	/////////////////////////////////////////////////////////////////////////////
	//right boundary
	double vel_u_right;
	double vel_v_right;
	//correction terms
	double correction_left;
	double correction_top_left;
	double correction_bottom_left;
  for(int j = 1; j < size_y - 1; j++)
  {
    int i = size_x - 1; //most right cell

		//compute missing values
		if (!useDirichletBc[3])
		{//Neumann -- outflow boundary
			//initialize as outflow
			double density = 1.0;
			//check if inflow - assume full pressure gradient is applied to inflow
			if (bcRight[0] > bcLeft[0])
			{
				density += 3 * (bcRight[0] - bcLeft[0]) / m_lattice_factor * size_x;
			}
			vel_u_right = (input[0](i, j) + input[2](i, j) + input[4](i, j) + 2 * (input[1](i, j) + input[5](i, j) + input[8](i, j))) - density;
			vel_v_right = 0;
		}
		else
		{//Dirichlet -- inflow or no-slip boundary
			vel_u_right = bcRight[0];
			vel_v_right = bcRight[1];
		}

		//set correction terms
		correction_left = - 2.0 / 3.0 * m_lattice_factor * vel_u_right;
		correction_top_left = - 0.5 * (input[2](i, j) - input[4](i, j)) + 0.5 * m_lattice_factor * vel_v_right - 1.0 / 6.0 * m_lattice_factor * vel_u_right;
		correction_bottom_left = 0.5 * (input[2](i, j) - input[4](i, j)) - 0.5 * m_lattice_factor * vel_v_right - 1.0 / 6.0 * m_lattice_factor * vel_u_right;;

    //half-way bounceback the values of the inner cell pointing out of the domain
    input[3](i, j) = input[1](i, j) + correction_left;
    input[6](i, j) = input[8](i, j) + correction_top_left;
    input[7](i, j) = input[5](i, j) + correction_bottom_left;

		//normal collide step
		collide(i, j, compute_rho(i,j, input), vel_u_right, vel_v_right, input);

		//modified stream step -> send values into domain
		//send input values to neighbor cells in output vector
		output[0](i, j) = input[0](i, j);
		output[2](i, j + 1) = input[2](i, j);
		output[3](i - 1, j) = input[3](i, j);
		output[4](i, j - 1) = input[4](i, j);
		output[6](i - 1, j + 1) = input[6](i, j);
		output[7](i - 1, j - 1) = input[7](i, j);
  }

	//////////////////////////////////////////////////////////////////////////////
	//left boundary
	double vel_u_left;
	double vel_v_left;
	//correction terms
	double correction_right;
	double correction_top_right;
	double correction_bottom_right;
  for(int j = 1; j < size_y - 1; j++)
  {
		int i = 0; //most left cell

		//compute missing values
		if (!useDirichletBc[2])
		{//Neumann -- inflow or outflow boundary
			//initialize as outflow
			double density = 1.0;
			//check if inflow - assume full pressure gradient is applied to inflow
			if (bcLeft[0] > bcRight[0])
			{
				density += 3 * (bcLeft[0] - bcRight[0]) / m_lattice_factor * size_x;
			}

			//compute velocity from pressure on boundary
			vel_u_left = density - (input[0](i, j) + input[2](i, j) + input[4](i, j) + 2 * (input[3](i, j) + input[6](i, j) + input[7](i, j)));
			vel_v_left = 0;
		}
		else
		{//Dirichlet -- inflow or no-slip boundary
			vel_u_left = bcLeft[0];
			vel_v_left = bcLeft[1];
		}

		//set correction terms
		correction_right =  2.0 / 3.0 * m_lattice_factor * vel_u_left;
		correction_top_right = - 0.5 * (input[2](i, j) - input[4](i, j)) + 0.5 * m_lattice_factor * vel_v_left + 1.0 / 6.0 * m_lattice_factor * vel_u_left;
		correction_bottom_right = 0.5 * (input[2](i, j) - input[4](i, j)) - 0.5 * m_lattice_factor * vel_v_left + 1.0 / 6.0 * m_lattice_factor * vel_u_left;

    //half-way bounceback the values of the inner cell pointing out of the domain
    input[1](i, j) = input[3](i, j) + correction_right;
    input[5](i, j) = input[7](i, j) + correction_top_right;
    input[8](i, j) = input[6](i, j) + correction_bottom_right;

		//normal collide step
		collide(i, j, compute_rho(i,j, input), vel_u_left, vel_v_left, input);

		//modified stream step -> send values into domain
		//send input values to neighbor cells in output vector
		output[0](i, j) = input[0](i, j);
		output[1](i + 1, j) = input[1](i, j);
		output[2](i, j + 1) = input[2](i, j);
		output[4](i, j - 1) = input[4](i, j);
		output[5](i + 1, j + 1) = input[5](i, j);
		output[8](i + 1, j - 1) = input[8](i, j);
	}

	//////////////////////////////////////////////////////////////////////////////
	//top boundary
	double vel_u_top;
	double vel_v_top;
	//correction terms
	double correction_bottom;
	double correction_left_bottom;
	double correction_right_bottom;
  for(int i = 1; i < size_x - 1; i++)
  {
    int j = size_y - 1; //most top cell

		//compute missing values
		if (!useDirichletBc[1])
		{//Neumann -- outflow boundary
			//initialize as outflow
			double density = 1.0;
			//check if inflow - assume full pressure gradient is applied to inflow
			if (bcTop[1] > bcBottom[1])
			{
				density += 3 * (bcTop[1] - bcBottom[1]) / m_lattice_factor * size_y;
			}

			//compute velocity from pressure on boundary
			vel_u_top = 0;
			vel_v_top = (input[0](i, j) + input[1](i, j) + input[3](i, j) + 2 * (input[2](i, j) + input[5](i, j) + input[6](i, j))) - density;
		}
		else
		{//Dirichlet -- inflow or no-slip boundary
			vel_u_top = bcTop[0];
			vel_v_top = bcTop[1];// * i * (size_x - 1 - i); //velocity profile doesnt work xD
		}

		//set correction terms
		correction_bottom =  - 2.0 / 3.0 * m_lattice_factor * vel_v_top;
		correction_left_bottom = 0.5 * (input[1](i, j) - input[3](i, j)) - 0.5 * m_lattice_factor * vel_u_top - 1.0 / 6.0 * m_lattice_factor * vel_v_top;
		correction_right_bottom = - 0.5 * (input[1](i, j) - input[3](i, j)) + 0.5 * m_lattice_factor * vel_u_top - 1.0 / 6.0 * m_lattice_factor * vel_v_top;

    //half-way bounceback the values of the inner cell pointing out of the domain
    input[4](i, j) = input[2](i, j) + correction_bottom;
    input[7](i, j) = input[5](i, j) + correction_left_bottom;
    input[8](i, j) = input[6](i, j) + correction_right_bottom;

		//normal collide step
		collide(i, j, compute_rho(i,j, input), vel_u_top, vel_v_top, input);

		//modified stream step -> send values into domain
		//send input values to neighbor cells in output vector
		output[0](i, j) = input[0](i, j);
		output[1](i + 1, j) = input[1](i, j);
		output[3](i - 1, j) = input[3](i, j);
		output[4](i, j - 1) = input[4](i, j);
		output[7](i - 1, j - 1) = input[7](i, j);
		output[8](i + 1, j - 1) = input[8](i, j);
  }

	//////////////////////////////////////////////////////////////////////////////
  //bottom boundary
	double vel_u_bottom;
	double vel_v_bottom;
	//correction terms
	double correction_top;
	double correction_right_top;
	double correction_left_top;
  for(int i = 1; i < size_x - 1; i++)
  {
    int j = 0; //most bottom cell

		//compute missing values
		if (!useDirichletBc[0])
		{//outflow boundary
			//initialize as outflow
			double density = 1.0;
			//check if inflow - assume full pressure gradient is applied to inflow
			if (bcBottom[1] > bcTop[1])
			{
				density += 3 * (bcBottom[1] - bcTop[1]) / m_lattice_factor * size_y;
			}

			//compute velocity from pressure on boundary
			vel_u_bottom = 0;
			vel_v_bottom = density - (input[0](i, j) + input[1](i, j) + input[3](i, j) + 2 * (input[4](i, j) + input[7](i, j) + input[8](i, j)));
		}
		else
		{//Dirichlet -- inflow or no-slip boundary
			vel_u_bottom = bcBottom[0];
			vel_v_bottom = bcBottom[1]; // * i * (size_x - 1 - i); //velocity profile doesnt work xD
		}

		//set correction terms
		correction_top = 2.0 / 3.0 * m_lattice_factor * vel_v_bottom;
		correction_right_top = - 0.5 * (input[1](i, j) - input[3](i, j)) + 0.5 * m_lattice_factor * vel_u_bottom + 1.0 / 6.0 * m_lattice_factor * vel_v_bottom;
		correction_left_top = 0.5 * (input[1](i, j) - input[3](i, j)) - 0.5 * m_lattice_factor * vel_u_bottom + 1.0 / 6.0 * m_lattice_factor * vel_v_bottom;

    //half-way bounceback the values of the inner cell pointing out of the domain
    input[2](i, j) = input[4](i, j) + correction_top;
    input[5](i, j) = input[7](i, j) + correction_right_top;
    input[6](i, j) = input[8](i, j) + correction_left_top;

		//normal collide step
		collide(i, j, compute_rho(i,j, input), vel_u_bottom , vel_v_bottom, input);

		//modified stream step -> send values into domain
		//send input values to neighbor cells in output vector
		output[0](i, j) = input[0](i, j);
		output[1](i + 1, j) = input[1](i, j);
		output[2](i, j + 1) = input[2](i, j);
		output[3](i - 1, j) = input[3](i, j);
		output[5](i + 1, j + 1) = input[5](i, j);
		output[6](i - 1, j + 1) = input[6](i, j);
  }
}
