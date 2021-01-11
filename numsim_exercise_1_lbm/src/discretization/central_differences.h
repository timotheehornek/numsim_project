#pragma once

#include "discretization.h"

class Central_differences :
	public Discretization
{
public:
	Central_differences(const std::array<int, 2>& nCells, const std::array<double, 2>& physicalSize, const double re, const std::array<double, 2>& g);

	//! compute the 1st derivative ∂ p / ∂x
	double compute_dp_dx(int i, int j) const override;

	//! compute the 2nd derivative ∂^2 u / ∂x^2
	double compute_du_dx2(int i, int j) const override;

	//! compute the 2nd derivative ∂^2 u / ∂y^2
	double compute_du_dy2(int i, int j) const override;

	//! compute the 1st derivative ∂ u^2 / ∂x
	double compute_du2_dx(int i, int j) const override;

	//! compute the 1st derivative ∂ (uv) / ∂y
	double compute_duv_dy(int i, int j) const override;


	//! compute the 1st derivative ∂ p / ∂y
	double compute_dp_dy(int i, int j) const override;

	//! compute the 2nd derivative ∂^2 v / ∂x^2
	double compute_dv_dx2(int i, int j) const override;

	//! compute the 2nd derivative ∂^2 v / ∂y^2
	double compute_dv_dy2(int i, int j) const override;

	//! compute the 1st derivative ∂ v^2 / ∂y
	double compute_dv2_dy(int i, int j) const override;

	//! compute the 1st derivative ∂ (uv) / ∂x
	double compute_duv_dx(int i, int j) const override;
};