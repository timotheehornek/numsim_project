#pragma once

#include "array2d/array2d.h"

#include <array>
#include <cassert>
#include <initializer_list>
#include <vector>

class Staggered_grid :
	public Array2D
{
private:
	enum Sides
	{
		SIDE_TOP,
		SIDE_LEFT,
		SIDE_RIGHT,
		SIDE_BOTTOM
	};
	enum Corners
	{
		CORNER_TOP_LEFT,
		CORNER_BOTTOM_RIGHT
	};
public:
	Staggered_grid(const std::array<int, 2>& size);
	//! constructor for initializer list
	Staggered_grid(const std::initializer_list<int>& l);

	//! return max index in x direction
	int x_max() const;

	//! return max index in y direction
	int y_max() const;

	//! return vector with inner boundary entries of current cell
	//! has distance 1 from array boundary
	//! (directions: 0=top, 1=left, 2=right, 3=bottom)
	std::vector<double> get_side(int side) const;

	//! write array boundary (without corners)
	//! (directions: 0=top, 1=left, 2=right, 3=bottom)
	void write_bound(int side, std::vector<double>& vec);

	//! return vector with every second inner boundary entry of current cell
	std::vector<double> get_side_odd(int side, int start) const;

	//! write array boundary (every second cell, without corners)
	void write_bound_odd(int side, std::vector<double> &vec, int start);

	
	//! return corners
	//! (corners: 0=top left, 1=bottom right)
	double get_corner(int corner) const;

	//! write corners
	//! (corners: 0=top left, 1=bottom right)
	void write_corner(int corner, const double val);

	//! interpolate values in grid (directions: 0=horizontal, 1=vertical, 2=horizontal+vertical)
	//! reference index is lower left or middle respectively for directions
	double interpolate(int i, int j, int direction) const;
};