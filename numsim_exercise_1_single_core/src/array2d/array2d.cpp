#include "array2d.h"



Array2D::Array2D(const std::array<int, 2>& size)
	: m_size{size}
{
	// allocate data, initialize to 0
	m_data.resize(m_size[0] * m_size[1], 0.0);
}

Array2D::Array2D(const std::initializer_list<int>& l)
	: m_size{ *l.begin(),*(l.begin()+1) }
{
	assert(l.size() == 2);

	// allocate data, initialize to 0
	m_data.resize(m_size[0] * m_size[1], 0.0);
}

//! get the size
std::array<int, 2> Array2D::size() const
{
	return m_size;
}

/*
double Array2D::interpolate(int i, int j, int direction) const
{
	assert(0 <= direction && direction <= 2);

	enum Directions{
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
*/
double Array2D::abs_max() const
{
	//! compute iterator to abs max of vector
	auto result_it = std::max_element(m_data.begin(), m_data.end(),
		[](double a, double b) {return (std::abs(a) < std::abs(b)); });
	return *result_it;
}

double& Array2D::operator()(int i, int j)
{
	const int index = j * m_size[0] + i;

	// assert that indices are in range
	assert(0 <= i && i < m_size[0]);
	assert(0 <= j && j < m_size[1]);
	assert(j * m_size[0] + i < (int)m_data.size());

	return m_data[index];
}

double Array2D::operator()(int i, int j) const
{
	const int index = j * m_size[0] + i;

	// assert that indices are in range
	assert(0 <= i && i < m_size[0]);
	assert(0 <= j && j < m_size[1]);
	assert(j * m_size[0] + i < (int)m_data.size());

	return m_data[index];
}

std::ostream& operator<<(std::ostream& out, const Array2D& array2d)
{
	out << array2d.m_size[0] << " X " << array2d.m_size[1] << " array:\n";
	for (int j{ array2d.m_size[1] - 1 }; j >= 0; --j)
	{
		for (int i{ 0 }; i < array2d.m_size[0]; ++i)
		{
			out << std::setprecision(3) << array2d(i, j) << '\t';
		}
		out << '\n';
	}
	return out;
}