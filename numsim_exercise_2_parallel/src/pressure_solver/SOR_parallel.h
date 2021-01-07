#pragma once

#include "pressure_solver.h"
#include <cmath>

class SOR_parallel
	: public Pressure_solver
{
private:
	const double m_w;

	const Discretization &m_disc;

	const std::array<int, 2> size;

	//! communicate with all neighbours / get boundary values of p
	void communicate_update_boundary_values(Staggered_grid &p, const bool black);

	//! iteration step on all red colored nodes (uneven indices)
	void step_black(Staggered_grid &p, const Staggered_grid &RHS);

	//! iteration step on all red colored nodes (even indices)
	void step_red(Staggered_grid &p, const Staggered_grid &RHS);

	std::vector<double> p_top_send_buffer{};
	std::vector<double> p_bottom_send_buffer{};
	std::vector<double> p_left_send_buffer{};
	std::vector<double> p_right_send_buffer{};

	std::vector<double> p_top_receive_buffer{};
	std::vector<double> p_bottom_receive_buffer{};
	std::vector<double> p_left_receive_buffer{};
	std::vector<double> p_right_receive_buffer{};

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
		TAG_RHS,
		TAG_SEND_TOP,
		TAG_SEND_LEFT,
		TAG_SEND_RIGHT,
		TAG_SEND_BOTTOM,

	};

public:
	// SOR_parallel(double dx, double dy, double eps, double max_it, double w);
	SOR_parallel(const Discretization &disc, Message_passer_2D &MP, double eps, int max_it, double w);

	// run one  Successive OverRelaxation iteration and update p
	void run_it_step(Staggered_grid &p, const Staggered_grid &RHS, bool black_red) override;
};
