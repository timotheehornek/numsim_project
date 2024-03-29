#pragma once

#include "discretization/staggered_grid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

class Lattice_boltzmann
{
protected:
  Staggered_grid m_u;
  Staggered_grid m_v;
  Staggered_grid m_p;

  Staggered_grid m_obstacle;

  double m_dt{ 0.0 }; //< run set_dt to initialize
  const double m_dx;
  const double m_dy;

  const double m_re;

  const double m_lattice_factor { 1.0 };
  const double m_time_factor { 0.01 };
  const double m_magic_factor;

  double m_one_over_tau{ 1.0 }; //< run set_tau to initialize

  const std::array<double, 2> m_g;

  const std::array<int, 2> m_nCells;

  //! compute hydrodynamic quantities
  double compute_rho(int i, int j, const std::array<Array2D, 9>& input) const;

  double compute_u(int i, int j, const std::array<Array2D, 9>& input) const;

  double compute_v(int i, int j, const std::array<Array2D, 9>& input) const;

  //! collide step at i,j
  void collide(int i, int j, double density, double vel_u, double vel_v, std::array<Array2D, 9>& input);

  //! bounceback values as object cell at i,j
  void bounceback(int i, int j, std::array<Array2D, 9>& input);

  //! stream i,j from input to output
  void stream(int i, int j, const std::array<Array2D, 9>& input, std::array<Array2D, 9>& output) const;

public:
  //! construct the object with given number of cells in x and y direction
	Lattice_boltzmann(const std::array<int, 2>& nCells, const std::array<double, 2>& physicalSize, const double re,
    const double magicFactor, const std::array<double, 2>& g);

	//! set dt without checking stability conditions
	void set_dt(double dt_max);

  //! set dt according to lattice factor
  void set_dt();

  //! set tau according to physical size and reynolds number
  void set_tau(double viscosity);

  //! set tau manually without the lattice factor
  void set_tau_manually(double tau);

 //! compute the viscosity of the fluid - relevant for setting tau
  double compute_viscosity(const std::array<double, 2>& size, const std::array<bool, 4>& useDirichletBc,
  		const std::array<double, 2>& bcBottom, const std::array<double, 2>& bcTop,
  		const std::array<double, 2>& bcLeft, const std::array<double, 2>& bcRight);

  //! set m_obstacle
  void set_obstacle_rect(const std::array<int, 4>& obstaclePos);

  void set_obstacle_round(const std::array<int, 2>& obstacleCenter, double radius);

  //! get max velocity set on boundary
  double get_max_vel(const std::array<bool, 4>& useDirichletBc,
    	const std::array<double, 2>& bcBottom, const std::array<double, 2>& bcTop,
    	const std::array<double, 2>& bcLeft, const std::array<double, 2>& bcRight) const;

  //! initialize f_vector in equilibrium state
  void initialize_f(double density, std::array<Array2D, 9>& input);

	//! getter
	double dt() const;
	double dx() const;
	double dy() const;
	const Staggered_grid& u() const;
	const Staggered_grid& v() const;
	const Staggered_grid& p() const;
	const std::array<int, 2>& nCells() const;

  const double u(int i, int j) const;
  const double v(int i, int j) const;
  const double p(int i, int j) const;

  //! run two lattice iteration with input and output on all inner nodes
  void two_lattice(std::array<Array2D, 9>& input, std::array<Array2D, 9>& output);

  //! boundary treatment for all boundary nodes
  void boundary_treatment(
  	std::array<Array2D, 9>& input, std::array<Array2D, 9>& output,
  	const std::array<bool, 4>& useDirichletBc,
  	const std::array<double, 2>& bcBottom, const std::array<double, 2>& bcTop,
  	const std::array<double, 2>& bcLeft, const std::array<double, 2>& bcRight);
};
