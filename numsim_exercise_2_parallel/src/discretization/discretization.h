#pragma once

#include "message_passing/message_passer.h"
#include "staggered_grid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

class Discretization
{
protected:
	//! DO NOT CHANGE ORDER OF INITIALIZATION!

	//! existing neighbors
	const std::array<bool, 4> m_ex_nbrs;
	const std::vector<int> m_nbrs;

	const std::array<int, 2> m_nCells;

	Staggered_grid m_u;
	Staggered_grid m_v;
	Staggered_grid m_p;
	Staggered_grid m_F;
	Staggered_grid m_G;
	Staggered_grid m_RHS;

	double m_dt{0.0}; //< run set_dt to initialize

	//! mesh width
	const double m_dx;
	const double m_dy;

	//! reynold number
	const double m_re;

	//! gravity
	const std::array<double, 2> m_g;

	Message_passer_2D &m_MP;

	enum Sides
	{
		SIDE_TOP,
		SIDE_LEFT,
		SIDE_RIGHT,
		SIDE_BOTTOM
	};
	enum Tags
	{
		TAG_U,
		TAG_V,
		TAG_F,
		TAG_G,
		TAG_P,
		TAG_RHS
	};
	enum Corners
	{
		CORNER_TOP_LEFT,
		CORNER_BOTTOM_RIGHT
	};

public:
	//! construct the object with given number of cells in x and y direction
	Discretization(const std::array<int, 2> &nCells_global, const std::array<double, 2> &physicalSize, const double re, const std::array<double, 2> &g, Message_passer_2D &MP);

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

	//! returns true if neighbor exists in direction
	//! (directions: 0=top, 1=left, 2=right, 3=bottom)
	bool ex_nbrs(int dir) const;

	//! returns ranks of neighbors (-1 if not existent)
	//! (directions: 0=top, 1=left, 2=right, 3=bottom)
	int nbrs(int dir) const;

	//! computes and updates boundary using boundary conditions of u,v around domain
	void communicate_update_bound_val_uv(
		const std::array<double, 2> &dirichletBcBottom, const std::array<double, 2> &dirichletBcTop,
		const std::array<double, 2> &dirichletBcLeft, const std::array<double, 2> &dirichletBcRight);

	//! compute and update F and G
	void compute_FG();

	//! communicate F and G required for RHS or update boundary values of F and G
	void communicate_update_bound_val_FG();

	//! compute and update RHS
	void compute_RHS();

	//! compute and update u and v
	void compute_uv();

	//! compute the 1st derivative ∂ p / ∂x
	virtual double compute_dp_dx(int i, int j) const = 0;

	//! compute the 2nd derivative ∂^2 u / ∂x^2
	virtual double compute_du_dx2(int i, int j) const = 0;

	//! compute the 2nd derivative ∂^2 u / ∂y^2
	virtual double compute_du_dy2(int i, int j) const = 0;

	//! compute the 1st derivative ∂ u^2 / ∂x
	virtual double compute_du2_dx(int i, int j) const = 0;

	//! compute the 1st derivative ∂ (uv) / ∂y
	virtual double compute_duv_dy(int i, int j) const = 0;

	//! compute the 1st derivative ∂ p / ∂y
	virtual double compute_dp_dy(int i, int j) const = 0;

	//! compute the 2nd derivative ∂^2 v / ∂x^2
	virtual double compute_dv_dx2(int i, int j) const = 0;

	//! compute the 2nd derivative ∂^2 v / ∂y^2
	virtual double compute_dv_dy2(int i, int j) const = 0;

	//! compute the 1st derivative ∂ v^2 / ∂y
	virtual double compute_dv2_dy(int i, int j) const = 0;

	//! compute the 1st derivative ∂ (uv) / ∂x
	virtual double compute_duv_dx(int i, int j) const = 0;
};