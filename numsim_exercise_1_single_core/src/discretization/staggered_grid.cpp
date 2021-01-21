#include "staggered_grid.h"

#include "central_differences.h"
#include "donor_cell.h"
#include "discretization.h"

Staggered_grid::Staggered_grid(const std::array<int, 2>& size)
	: Array2D{ size }, m_x_max{ size[0] }, m_y_max{ size[1] }
{}

Staggered_grid::Staggered_grid(const std::initializer_list<int>& l)
	: Array2D{ l }, m_x_max{ *l.begin() }, m_y_max{ *(l.begin()+1) }
{}


//! copy constructor
Staggered_grid::Staggered_grid(const Staggered_grid &copy)
	: Array2D{ copy.m_size }, m_x_max{ copy.m_size[0] }, m_y_max{ copy.m_size[1] }
{}

//! overloaded assignment
Staggered_grid& Staggered_grid::operator=(const Staggered_grid& copy)
{
	// self-assignment check
	if (this == &copy)
		return *this;
	
	assert(m_size[0]==copy.m_size[0]&&m_size[1]==copy.m_size[1]);
	
	m_data=copy.m_data;
	return *this;
}

double Staggered_grid::x_max() const
{
	return m_size[0];
}

double Staggered_grid::y_max() const
{
	return m_size[1];
}

double Staggered_grid::interpolate(int i, int j, int direction) const
{
	assert(0 <= direction && direction <= 2);

	enum Directions {
		HORIZONTAL,
		VERTICAL,
		HORIZONTAL_VERTICAL,
	};

	switch (direction) {
	case HORIZONTAL:
		return ((*this)(i, j) + (*this)(i + 1, j)) / 2.0; break;
	case VERTICAL:
		return ((*this)(i, j) + (*this)(i, j + 1)) / 2.0; break;
	case HORIZONTAL_VERTICAL:
		return ((*this)(i, j) + (*this)(i + 1, j) + (*this)(i, j + 1) + (*this)(i + 1, j + 1)) / 4.0; break;
	default:
		return 0.0;
	}
}