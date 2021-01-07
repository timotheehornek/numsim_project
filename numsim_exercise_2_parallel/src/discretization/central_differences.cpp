#include "central_differences.h"

Central_differences::Central_differences(const std::array<int, 2> &nCells, const std::array<double, 2> &physicalSize, const double re, const std::array<double, 2> &g, Message_passer_2D &MP)
	: Discretization(nCells, physicalSize, re, g, MP){};

//! compute the 1st derivative ∂ p / ∂x
double Central_differences::compute_dp_dx(int i, int j) const
{
	return (m_p(i + 1, j) - m_p(i, j)) / m_dx;
}

//! compute the 2nd derivative ∂^2 u / ∂x^2
double Central_differences::compute_du_dx2(int i, int j) const
{
	//! compute central differences u_x_i_plushalf_j and u_x_i_minushalf_j
	//! and use these to compute central difference quotient u_xx_i_j
	return (m_u(i + 1, j) - 2. * m_u(i, j) + m_u(i - 1, j)) / (m_dx * m_dx);
}

//! compute the 2nd derivative ∂^2 u / ∂y^2
double Central_differences::compute_du_dy2(int i, int j) const
{
	//! compute central differences u_y_i_jplushalf and u_y_i_jminushalf
	//! and use these to compute central difference quotient u_yy_i_j
	return (m_u(i, j + 1) - 2. * m_u(i, j) + m_u(i, j - 1)) / (m_dy * m_dy);
}

//! compute the 1st derivative ∂ u^2 / ∂x
double Central_differences::compute_du2_dx(int i, int j) const
{
	//! compute u_iplushalf_j and u_iminushalf_j by using linear interpolation between their neighbours,
	//! square them and compute central differences
	const double u_iplushalf_j = 1. / 2. * (m_u(i + 1, j) + m_u(i, j));
	const double u_iminushalf_j = 1. / 2. * (m_u(i, j) + m_u(i - 1, j));
	return ((u_iplushalf_j * u_iplushalf_j) - (u_iminushalf_j * u_iminushalf_j)) / m_dx;
}

//! compute the 1st derivative ∂ (uv) / ∂y
double Central_differences::compute_duv_dy(int i, int j) const
{
	//! compute u_i_jplushalf, v_i_jplushalf, u_i_jminushalf, v_i_jminushalf by using linear interpolation //! between their neighbours
	//! and compute central differences
	const double uv_i_jplushalf = ((m_u(i, j + 1) + m_u(i, j)) * (m_v(i, j) + m_v(i + 1, j))) / 4.;			 // uv at top right corner of cell
	const double uv_i_jminushalf = ((m_u(i, j) + m_u(i, j - 1)) * (m_v(i, j - 1) + m_v(i + 1, j - 1))) / 4.; // uv at bottom right corner of cell
	return (uv_i_jplushalf - uv_i_jminushalf) / m_dy;
}

//! compute the 1st derivative ∂ p / ∂y
double Central_differences::compute_dp_dy(int i, int j) const
{
	return (m_p(i, j + 1) - m_p(i, j)) / m_dy;
}

//! compute the 2nd derivative ∂^2 v / ∂x^2
double Central_differences::compute_dv_dx2(int i, int j) const
{
	//! compute central differences v_x_i_plushalf_j and v_x_i_minushalf_j
	//! and use these to compute central difference qoutient v_xx_i_j
	return (m_v(i + 1, j) - 2. * m_v(i, j) + m_v(i - 1, j)) / (m_dx * m_dx);
}

//! compute the 2nd derivative ∂^2 v / ∂y^2
double Central_differences::compute_dv_dy2(int i, int j) const
{
	//! compute central differences v_y_i_jplushalf and v_y_i_jminushalf
	//! and use these to compute central difference quotient v_yy_i_j
	return (m_v(i, j + 1) - 2. * m_v(i, j) + m_v(i, j - 1)) / (m_dy * m_dy);
}

//! compute the 1st derivative ∂ v^2 / ∂y
double Central_differences::compute_dv2_dy(int i, int j) const
{
	//! compute v_i_jplushalf and v_i_jminushalf by using linear interpolation between their neighbours,
	//! square them and compute central differences
	const double u_i_jplushalf = 1. / 2. * (m_v(i, j + 1) + m_v(i, j));
	const double u_i_jminushalf = 1. / 2. * (m_v(i, j) + m_v(i, j - 1));
	return ((u_i_jplushalf * u_i_jplushalf) - (u_i_jminushalf * u_i_jminushalf)) / m_dy;
}

//! compute the 1st derivative ∂ (uv) / ∂x
double Central_differences::compute_duv_dx(int i, int j) const
{
	//! compute u_iplushalf_j, v_iplushalf_j, u_iminushalf_j, v_iminushalf_j by using linear interpolation //! between their neighbours
	//! and compute central differences
	const double uv_iplushalf_j = ((m_u(i, j + 1) + m_u(i, j)) * (m_v(i, j) + m_v(i + 1, j))) / 4.;			 // uv at top right corner of cell
	const double uv_iminushalf_j = ((m_u(i - 1, j + 1) + m_u(i - 1, j)) * (m_v(i, j) + m_v(i - 1, j))) / 4.; // uv at top left corner of cell
	return (uv_iplushalf_j - uv_iminushalf_j) / m_dx;
}