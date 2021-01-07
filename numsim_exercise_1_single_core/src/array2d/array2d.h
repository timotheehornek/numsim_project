#pragma once

#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <string>
#include <vector>

/** This class represents a 2D array of double values.
 *  Internally they are stored consecutively in memory.
 *  The entries can be accessed by two indices i,j.
 */
class Array2D
{
protected:
	std::vector<double> m_data;  //< storage array values, in row-major order
	const std::array<int, 2> m_size;    //< width, height of the domain

public:
	//! constructor
	Array2D(const std::array<int, 2>& size);

	//! constructor for initializer list
	Array2D(const std::initializer_list<int>& l);

	//! get the size
	std::array<int, 2> size() const;

	//! interpolate values in array (directions: 0=horizontal, 1=vertical, 2=horizontal+vertical)
	//double interpolate(int i, int j, int direction) const;

	//! get largest absolute entry
	double abs_max() const;

	//! access the value at coordinate (i,j), declared not const, i.e. the value can be changed
	double& operator()(int i, int j);

	//! get the value at coordinate (i,j), declared const, i.e. it is not possible to change the value
	double operator()(int i, int j) const;

	//! outputstream of 2D array
	friend std::ostream& operator<< (std::ostream& out, const Array2D& array2d);
};