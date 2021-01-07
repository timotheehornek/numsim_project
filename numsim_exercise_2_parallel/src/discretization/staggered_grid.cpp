#include "staggered_grid.h"

Staggered_grid::Staggered_grid(const std::array<int, 2> &size)
	: Array2D{size}{}

Staggered_grid::Staggered_grid(const std::initializer_list<int> &l)
	: Array2D{l}{}

int Staggered_grid::x_max() const
{
	return m_size[0];
}

int Staggered_grid::y_max() const
{
	return m_size[1];
}

std::vector<double> Staggered_grid::get_side(int side) const
{
	assert(side >= 0 && side < 4);

	std::vector<double> vec{};

	switch (side)
	{
	case SIDE_TOP:
		vec.reserve(m_size[0] - 2);
		for (int i{1}; i < m_size[0] - 1; ++i)
		{
			vec.push_back((*this)(i, m_size[1] - 2));
		}
		break;
	case SIDE_LEFT:
		vec.reserve(m_size[1] - 2);
		for (int j{1}; j < m_size[1] - 1; ++j)
		{
			vec.push_back((*this)(1, j));
		}
		break;
	case SIDE_RIGHT:
		vec.reserve(m_size[1] - 2);
		for (int j{1}; j < m_size[1] - 1; ++j)
		{
			vec.push_back((*this)(m_size[0] - 2, j));
		}
		break;
	case SIDE_BOTTOM:
		vec.reserve(m_size[0] - 2);
		for (int i{1}; i < m_size[0] - 1; ++i)
		{
			vec.push_back((*this)(i, 1));
		}
		break;
	}
	return vec;
}

void Staggered_grid::write_bound(int side, std::vector<double> &vec)
{
	assert(side >= 0 && side < 4);

	switch (side)
	{
	case SIDE_TOP:
		assert(vec.size() == m_size[0] - 2);
		for (int i{0}; i < m_size[0] - 2; ++i)
		{
			(*this)(i + 1, m_size[1] - 1) = vec[i];
		}
		break;
	case SIDE_LEFT:
		assert(vec.size() == m_size[1] - 2);
		for (int j{0}; j < m_size[1] - 2; ++j)
		{
			(*this)(0, j + 1) = vec[j];
		}
		break;
	case SIDE_RIGHT:
		assert(vec.size() == m_size[1] - 2);
		for (int j{0}; j < m_size[1] - 2; ++j)
		{
			(*this)(m_size[0] - 1, j + 1) = vec[j];
		}
		break;
	case SIDE_BOTTOM:
		assert(vec.size() == m_size[0] - 2);
		for (int i{}; i < m_size[0] - 2; ++i)
		{
			(*this)(i + 1, 0) = vec[i];
		}
		break;
	}
}

std::vector<double> Staggered_grid::get_side_odd(int side, int start) const
{
	assert(side >= 0 && side < 4);

	std::vector<double> vec{};

	switch (side)
	{
	case SIDE_TOP:
		vec.reserve(m_size[0] /2);
		for (int i{1+start}; i < m_size[0] - 1; i+=2)
		{
			vec.push_back((*this)(i, m_size[1] - 2));
		}
		break;
	case SIDE_LEFT:
		vec.reserve(m_size[1] /2);
		for (int j{1+start}; j < m_size[1] - 1; j+=2)
		{
			vec.push_back((*this)(1, j));
		}
		break;
	case SIDE_RIGHT:
		vec.reserve(m_size[1] /2);
		for (int j{1+start}; j < m_size[1] - 1; j+=2)
		{
			vec.push_back((*this)(m_size[0] - 2, j));
		}
		break;
	case SIDE_BOTTOM:
		vec.reserve(m_size[0] /2);
		for (int i{1+start}; i < m_size[0] - 1; i+=2)
		{
			vec.push_back((*this)(i, 1));
		}
		break;
	}
	return vec;
}

void Staggered_grid::write_bound_odd(int side, std::vector<double> &vec, int start)
{
	assert(side >= 0 && side < 4);
	//! define vector index (different from write index)
	int vec_idx{0};
	switch (side)
	{
	case SIDE_TOP:
		assert(vec.size() == m_size[0] /2);
		for (int i{start}; i < m_size[0] - 2; i+=2, ++vec_idx)
		{
			(*this)(i + 1, m_size[1] - 1) = vec[vec_idx];
		}
		break;
	case SIDE_LEFT:
		assert(vec.size() == m_size[1] /2);
		for (int j{start}; j < m_size[1] - 2; j+=2, ++vec_idx)
		{
			(*this)(0, j + 1) = vec[vec_idx];
		}
		break;
	case SIDE_RIGHT:
		assert(vec.size() == m_size[1] /2);
		for (int j{start}; j < m_size[1] - 2; j+=2, ++vec_idx)
		{
			(*this)(m_size[0] - 1, j + 1) = vec[vec_idx];
		}
		break;
	case SIDE_BOTTOM:
		assert(vec.size() == m_size[0] /2);
		for (int i{start}; i < m_size[0] - 2; i+=2, ++vec_idx)
		{
			(*this)(i + 1, 0) = vec[vec_idx];
		}
		break;
	}
}

double Staggered_grid::get_corner(int corner) const
{
	assert(0 <= corner && corner < 2);

	double val{};
	switch (corner)
	{
	case CORNER_TOP_LEFT:
		val = (*this)(1, m_size[1] - 2);
		break;
	case CORNER_BOTTOM_RIGHT:
		val = (*this)(m_size[0] - 2, 1);
		break;
	}
	return val;
}

void Staggered_grid::write_corner(int corner, const double val)
{
	assert(0 <= corner && corner < 2);

	switch (corner)
	{
	case CORNER_TOP_LEFT:
		(*this)(0, m_size[1] - 1) = val;
		break;
	case CORNER_BOTTOM_RIGHT:
		(*this)(m_size[0] - 1, 0) = val;
		break;
	}
}

double Staggered_grid::interpolate(int i, int j, int direction) const
{
	assert(0 <= direction && direction <= 2);

	enum Directions
	{
		HORIZONTAL,
		VERTICAL,
		HORIZONTAL_VERTICAL,
	};

	switch (direction)
	{
	case HORIZONTAL:
		return ((*this)(i, j) + (*this)(i + 1, j)) / 2.0;
		break;
	case VERTICAL:
		return ((*this)(i, j) + (*this)(i, j + 1)) / 2.0;
		break;
	case HORIZONTAL_VERTICAL:
		return ((*this)(i, j) + (*this)(i + 1, j) + (*this)(i, j + 1) + (*this)(i + 1, j + 1)) / 4.0;
		break;
	default:
		return 0.0;
	}
}