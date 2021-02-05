#pragma once

#include "pressure_solver.h"
#include "discretization/staggered_grid.h"
#include <cmath>

class CG
	: public Pressure_solver
{
private:
	// compute all variables needed for CG algorithm
	// A == System-matrix
	// p == search_directions; also called "directions" to avoid confusion with pressure
	// r == residual
	// alpha_k, beta_k == auxiliary variables needed for computation in the k-th step

	//! pressure array size
	std::array<int, 2> m_size;

	//! staggered_grids for auxiliary variables
	Staggered_grid m_residuals;
	Staggered_grid m_directions;
	Staggered_grid m_Ap;


	// Matrix-vector-product
	void compute_Ap(Discretization &discr);
	// squared L2-norm of vector
	double compute_rk_transp_rk(const Staggered_grid& residuals);
	// squared L2-norm of p^T*A*p
	double compute_pk_transp_Apk();
	// set new pressure values
	void update_p(double alpha_k, Discretization &discr);
	// compute new residual
	void update_residual(double alpha_k);
	// compute new search directions
	void update_directions(double beta_k);

	//! compute residual using p and RHS for each gridpoint; store value in "residuals"
	void compute_residual(Discretization &discr);

	// run one  conjugated gradient iteration and update p
	void run_it_step(Discretization &discr);

	enum SIDES
	{
		BOTTOM,
		TOP,
		LEFT,
		RIGHT
	};

//protected:
	

	//store all values of send in recv
	//void copy_array(const Staggered_grid &send, Staggered_grid &recv, Discretization &discr) override;

public:
	// initialize CG-algorithm
	CG(double eps, int max_it, std::array<int, 2> nCells);

	// implementation of CG pressure solver
	void solver(Discretization &discr,
				const std::array<double, 2> &bcBottom, const std::array<double, 2> &bcTop,
				const std::array<double, 2> &bcLeft, const std::array<double, 2> &bcRight,
				const std::array<bool, 4> &useDirichletBc) override;
				
};
