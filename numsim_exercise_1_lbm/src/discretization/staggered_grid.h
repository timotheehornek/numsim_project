#pragma once

#include "array2d/array2d.h"

#include <array>
#include <cassert>
#include <initializer_list>

class Staggered_grid :
	public Array2D
{
private:
	const int m_x_max;
	const int m_y_max;

public:
	Staggered_grid(const std::array<int, 2>& size);
	//! constructor for initializer list
	Staggered_grid(const std::initializer_list<int>& l);

	//! return max index in x direction
	double x_max() const;

	//! return max index in x direction
	double y_max() const;

	//! interpolate values in grid (directions: 0=horizontal, 1=vertical, 2=horizontal+vertical)
	//! reference index is lower left or middle respectively for directions
	double interpolate(int i, int j, int direction) const;
};