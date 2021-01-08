#include "SOR.h"

SOR::SOR(double dx, double dy, double eps, double max_it,double w)
	: Pressure_solver(dx, dy, eps, max_it), m_w{w}{}

//void SOR::run_it_step(Array2D& p, const Array2D& RHS, const std::array<int, 2> size) const
void SOR::run_it_step(Discretization& discr) const override;
{
	for (int j{ 1 }; j < discr.p().size()[1] - 1; ++j)
	{
		for (int i{ 1 }; i < discr.p().size()[0] - 1; ++i)
		{
			//! single node iteration step
			discr.p_ref(i, j) =
				(1-m_w)*p(i, j)
				+ m_w																	   //< factor omega
				* std::pow(m_dx * m_dy, 2) / (2 * (std::pow(m_dx, 2) + std::pow(m_dy, 2))) //< prefactor
				* ((discr.p(i - 1, j) + discr.p(i + 1, j)) / std::pow(m_dx, 2)						   //< horizontal term
					+ (discr.p(i, j - 1) + discr.p(i, j + 1)) / std::pow(m_dy, 2)					   //< vertical term
					- discr.RHS(i, j));															//< right hand side
																		       
		}
		//! update boundary values on left and right
		discr.p_ref(0, j) = p(1, j);
		discr.p_ref(discr.p().size()[0] - 1, j) = p(discr.p().size()[0] - 2, j);
	}

	//! update boundary values on bottom and top
	for (int i{ 0 }; i < discr.p().size()[0] ; ++i)
	{
		discr.p_ref(i, 0) = discr.p(i, 1);
		discr.p_ref(i, discr.p().size()[1] - 1) = discr.p(i, discr.p().size()[1] - 2);
	}
}