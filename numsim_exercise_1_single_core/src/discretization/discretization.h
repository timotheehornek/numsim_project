#pragma once

#include "staggered_grid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

class Discretization
{
protected:
	Staggered_grid m_u;
	Staggered_grid m_v;
	Staggered_grid m_p;
	Staggered_grid m_F;
	Staggered_grid m_G;
	Staggered_grid m_RHS;

	double m_dt{0.0}; //< run set_dt to initialize
	const double m_dx;
	const double m_dy;

	const double m_re;

	const std::array<double, 2> m_g;

	const std::array<int, 2> m_nCells;

public:
	//! construct the object with given number of cells in x and y direction
	Discretization(const std::array<int, 2> &nCells, const std::array<double, 2> &physicalSize, const double re, const std::array<double, 2> &g);

	//! set dt checking stability conditions
	void set_dt(double dt_max, double sec_factor);
	//! set dt without checking stability conditions
	void set_dt(double dt_max);

	//! getter
	double dt() const;
	double dx() const;
	double dy() const;
	const Staggered_grid &u() const;
	const Staggered_grid &v() const;
	const Staggered_grid &p() const;
	const Staggered_grid &F() const;
	const Staggered_grid &G() const;
	const Staggered_grid &RHS() const;
	Staggered_grid &p_ref(); //< return p by non-const reference
	const std::array<int, 2> &nCells() const;

	//! sets up boundary using boundary conditions of u,v around domain
	void setup_bound_val_uv(
		const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
		const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
		const std::array<bool, 4> &useDirichletBc);

	//! computes and updates boundary using boundary conditions of u,v around domain
	void update_bound_val_uv(
		const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
		const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
		const std::array<bool, 4> &useDirichletBc);

	//! compute and update boundary values of F and G
	void compute_bound_val_FG();

	// compute and update F and G
	void compute_FG();

	//compute and update RHS
	void compute_RHS();

	//compute and update u and v
	void compute_uv();

	// compute the 1st derivative ∂ p / ∂x
	virtual double compute_dp_dx(int i, int j) const = 0;

	// compute the 2nd derivative ∂^2 u / ∂x^2
	virtual double compute_du_dx2(int i, int j) const = 0;

	// compute the 2nd derivative ∂^2 u / ∂y^2
	virtual double compute_du_dy2(int i, int j) const = 0;

	// compute the 1st derivative ∂ u^2 / ∂x
	virtual double compute_du2_dx(int i, int j) const = 0;

	// compute the 1st derivative ∂ (uv) / ∂y
	virtual double compute_duv_dy(int i, int j) const = 0;

	// compute the 1st derivative ∂ p / ∂y
	virtual double compute_dp_dy(int i, int j) const = 0;

	// compute the 2nd derivative ∂^2 v / ∂x^2
	virtual double compute_dv_dx2(int i, int j) const = 0;

	// compute the 2nd derivative ∂^2 v / ∂y^2
	virtual double compute_dv_dy2(int i, int j) const = 0;

	// compute the 1st derivative ∂ v^2 / ∂y
	virtual double compute_dv2_dy(int i, int j) const = 0;

	// compute the 1st derivative ∂ (uv) / ∂x
	virtual double compute_duv_dx(int i, int j) const = 0;
};