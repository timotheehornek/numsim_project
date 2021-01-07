#include "donor_cell.h"

Donor_cell::Donor_cell(
	const std::array<int, 2>& nCells, 
	const std::array<double, 2>& physicalSize, 
	const double re, const std::array<double, 2>& g, 
	Message_passer_2D& MP, 
	const double alpha)
	: Discretization(nCells, physicalSize, re, g, MP), m_alpha{alpha}
	{}

//! compute the 1st derivative ∂ p / ∂x
double Donor_cell::compute_dp_dx(int i, int j) const
{
	return (m_p(i + 1, j) - m_p(i, j)) / m_dx;
}

//! compute the 2nd derivative ∂^2 u / ∂x^2
double Donor_cell::compute_du_dx2(int i, int j) const
{
	//! compute central differences u_x_i_plushalf_j and u_x_i_minushalf_j
	//! and use these to compute central difference quotient u_xx_i_j
	return (m_u(i + 1, j) - 2. * m_u(i, j) + m_u(i - 1, j)) / (m_dx * m_dx);
}

//! compute the 2nd derivative ∂^2 u / ∂y^2
double Donor_cell::compute_du_dy2(int i, int j) const
{
	//! compute central differences u_y_i_jplushalf and u_y_i_jminushalf
    //! and use these to compute central difference quotient u_yy_i_j
	return (m_u(i, j + 1) - 2. * m_u(i, j) + m_u(i, j - 1)) / (m_dy * m_dy);
}

//! compute the 1st derivative ∂ u^2 / ∂x
double Donor_cell::compute_du2_dx(int i, int j) const
{
	//! compute u_iplushalf_j and u_iminushalf_j by using linear interpolation between their neighbours,
	//! square them and compute central differences,
	//! then add the donor cell upwind scheme
	const double u_iplushalf_j = 1. / 2. * (m_u(i + 1, j) + m_u(i, j));
	const double u_iminushalf_j = 1. / 2. * (m_u(i, j) + m_u(i - 1, j));
	const double central_difference = (u_iplushalf_j * u_iplushalf_j - u_iminushalf_j * u_iminushalf_j) / m_dx;

	const double donor_u_iplushalf_j = 1. / 2. * (m_u(i, j) - m_u(i + 1, j));
	const double donor_u_iminushalf_j = 1. / 2. * (m_u(i - 1, j) - m_u(i, j)); //const double donor_u_iminushalf_j = 1. / 2. * (m_u(i - 1, j) + m_u(i, j));
	const double donor_cell_modification = (std::abs(u_iplushalf_j) * donor_u_iplushalf_j - std::abs(u_iminushalf_j) * donor_u_iminushalf_j) / m_dx;

	return central_difference + m_alpha * donor_cell_modification;
}

//! compute the 1st derivative ∂ (uv) / ∂y
double Donor_cell::compute_duv_dy(int i, int j) const
{
	//! compute u_i_jplushalf, v_iplushalf_j, u_i_jminushalf, v_iplushalf_jminus by using linear interpolation between their neighbours
	//! and compute central differences
	const double u_i_jplushalf = (m_u(i, j + 1) + m_u(i, j)) / 2.;       // u at top right corner of cell
	const double v_iplushalf_j = (m_v(i, j) + m_v(i + 1, j)) / 2.;     // v at top right corner of cell
	const double u_i_jminushalf = (m_u(i, j) + m_u(i, j - 1)) / 2.;     // u at bottom right corner of cell
	const double v_iplushalf_jminus = (m_v(i, j - 1) + m_v(i + 1, j - 1)) / 2.;    // v at bottom right corner of cell

	const double central_difference = (v_iplushalf_j * u_i_jplushalf - v_iplushalf_jminus * u_i_jminushalf) / m_dy;

	//! then compute the upwind scheme for the donor cell modification
	const double donor_u_i_jplushalf = (m_u(i, j) - m_u(i, j + 1)) / 2.; // for looking upward flow direction at upper corner of cell
	const double donor_u_i_jminushalf = (m_u(i, j - 1) - m_u(i, j)) / 2.; // for looking upward flow direction at lower corner of cell

	const double donor_cell_modification = (std::abs(v_iplushalf_j) * donor_u_i_jplushalf - std::abs(v_iplushalf_jminus) * donor_u_i_jminushalf) / m_dy;

	//! compute the weighted output for the donor cell
	return central_difference + m_alpha * donor_cell_modification;
}


//! compute the 1st derivative ∂ p / ∂y
double Donor_cell::compute_dp_dy(int i, int j) const
{
	return (m_p(i, j + 1) - m_p(i, j)) / m_dy;
}

//! compute the 2nd derivative ∂^2 v / ∂x^2
double Donor_cell::compute_dv_dx2(int i, int j) const
{
	//! compute central differences v_x_i_plushalf_j and v_x_i_minushalf_j
	//! and use these to compute central difference qoutient v_xx_i_j
	return (m_v(i + 1, j) - 2. * m_v(i, j) + m_v(i - 1, j)) / (m_dx * m_dx);
}

//! compute the 2nd derivative ∂^2 v / ∂y^2
double Donor_cell::compute_dv_dy2(int i, int j) const
{
	//! compute central differences v_y_i_jplushalf and v_y_i_jminushalf
    //! and use these to compute central difference quotient v_yy_i_j
	return (m_v(i, j + 1) - 2. * m_v(i, j) + m_v(i, j - 1)) / (m_dy * m_dy);
}

//! compute the 1st derivative ∂ v^2 / ∂y
double Donor_cell::compute_dv2_dy(int i, int j) const
{
	//! compute u_iplushalf_j and u_iminushalf_j by using linear interpolation between their neighbours,
   //! square them and compute central differences,
   //! then add the donor cell upwind scheme
	const double v_i_jplushalf = 1. / 2. * (m_v(i, j + 1) + m_v(i, j));
	const double v_i_jminushalf = 1. / 2. * (m_v(i, j) + m_v(i, j - 1));
	const double central_difference = (v_i_jplushalf * v_i_jplushalf - v_i_jminushalf * v_i_jminushalf) / m_dy;

	const double donor_v_i_jplushalf = 1. / 2. * (m_v(i, j) - m_v(i, j + 1));
	const double donor_v_i_jminushalf = 1. / 2. * (m_v(i, j - 1) - m_v(i, j));
	const double donor_cell_modification = (std::abs(v_i_jplushalf) * donor_v_i_jplushalf - std::abs(v_i_jminushalf) * donor_v_i_jminushalf) / m_dy;

	return central_difference + m_alpha * donor_cell_modification;
}

//! compute the 1st derivative ∂ (uv) / ∂x
double Donor_cell::compute_duv_dx(int i, int j) const
{
	//! compute u_i_jplushalf, v_iplushalf_j, u_i_jminushalf, v_iplushalf_jminus by using linear interpolation between their neighbours
	//! and compute central differences
	const double u_i_jplushalf = (m_u(i, j + 1) + m_u(i, j)) / 2.;       // u at top right corner of cell
	const double v_iplushalf_j = (m_v(i, j) + m_v(i + 1, j)) / 2.;     // v at top right corner of cell
	const double u_iminus_jplushalf = (m_u(i - 1, j) + m_u(i - 1, j + 1)) / 2.;     // u at top left corner of cell
	const double v_iminushalf_j = (m_v(i - 1, j) + m_v(i, j)) / 2.;    // v at top left corner of cell

	const double central_difference = (u_i_jplushalf * v_iplushalf_j - u_iminus_jplushalf * v_iminushalf_j) / m_dx;

	//! then compute the upwind scheme for the donor cell modification
	const double donor_v_iplushalf_j = (m_v(i, j) - m_v(i + 1, j)) / 2.; // for looking upward flow direction at upper corner of cell
	const double donor_v_iminushalf_j = (m_v(i - 1, j) - m_v(i, j)) / 2.; // for looking upward flow direction at lower corner of cell

	const double donor_cell_modification = (std::abs(u_i_jplushalf) * donor_v_iplushalf_j - std::abs(u_iminus_jplushalf) * donor_v_iminushalf_j) / m_dx;

	//! compute the weighted output for the donor cell
	return central_difference + m_alpha * donor_cell_modification;
}
