// Author: Damian Marek

// Copyright 2021 Shashwat Sharma and Piero Triverio

// This file is part of Strata.

// Strata is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Strata is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Strata.  If not, see <https://www.gnu.org/licenses/>.

/************************ interpolation_table.hpp ************************

 * Interpolation table for 3D datasets using piecewise polynomials.
 *
 * Author: Damian Marek
 * Created on: April 09, 2022

 ***************************************************************/

#ifndef STRATA_INTERPOLATION_TABLE_HPP
#define STRATA_INTERPOLATION_TABLE_HPP

#include <cassert>
#include <functional>

using Int3D = std::array<int, 3>;
using Double3D = std::array<double, 3>;

inline Int3D operator-(const Int3D &lhs, const Int3D &rhs)
{
	Int3D diff;
	for (int i = 0; i < 3; i++)
		diff[i] = lhs[i] - rhs[i];
	return diff;
}

inline Int3D operator+(const Int3D &lhs, const Int3D &rhs)
{
	Int3D sum;
	for (int i = 0; i < 3; i++)
		sum[i] = lhs[i] + rhs[i];
	return sum;
}

inline Int3D operator+(const Int3D &lhs, const int rhs)
{
	Int3D sum;
	for (int i = 0; i < 3; i++)
		sum[i] = lhs[i] + rhs;
	return sum;
}
template <size_t N>
inline std::array<std::complex<double>, N> operator*(const std::array<std::complex<double>, N> &lhs,
													 const double rhs)
{
	std::array<std::complex<double>, N> result;
	for (int i = 0; i < N; i++)
		result[i] = lhs[i] * rhs;
	return result;
}
template <size_t N>
inline std::array<std::complex<double>, N> &
operator+=(std::array<std::complex<double>, N> &lhs, const std::array<std::complex<double>, N> &rhs)
{
	for (int i = 0; i < N; i++)
		lhs[i] += rhs[i];
	return lhs;
}

template <typename Iterator, typename T>
ptrdiff_t find_closest_point(T item, Iterator begin, Iterator end)
{
	assert(std::is_sorted(begin, end));
	assert(*begin <= item);
	assert(item <= *(end - 1));
	// check edge cases
	if (*begin == item)
		return 0;
	else if (item == *(end - 1))
		return end - 1 - begin;
	auto lb = std::lower_bound(begin, end, item);
	assert(lb != end);
	if (*lb == item)
	{
		// test point is coincident with a data point
		return lb - begin;
	}
	else // test < *lb
	{
		T other_b = *(lb - 1);
		return std::abs(other_b - item) < std::abs(*lb - item) ? lb - begin - 1 : lb - begin;
	}
}
/// Determine the interval where an item belongs.
/// Each interval is defined by a range [a, b) where "b" is not included
/// A collection of ranges are defined as {a, b, c, ...}, where the first interval is [a, b), and
/// the second is defined as [b, c), etc. The final interval will include its end point in the
/// interval.
template <typename Iterator, typename T>
ptrdiff_t find_interval(T item, Iterator begin, Iterator end)
{
	assert(std::is_sorted(begin, end));
	assert(*begin <= item);
	assert(item <= *(end - 1));
	// check edge cases
	if (*begin == item)
		return 0;
	else if (item == *(end - 1))
		return end - 1 - begin;

	auto ub = std::upper_bound(begin, end, item);
	assert(ub != end);
	ptrdiff_t interval = (ub - begin - 1);
	return interval;
}

/// Piecewise polynomial interpolation of 3D data
template <int N> class InterpolationTable
{
  public:
	/// Sets the grid and the associated data. Data is stored row-major, so the "z" coordinate
	/// changes the quickest.
	InterpolationTable(std::vector<double> _x, std::vector<double> _y, std::vector<double> _z,
					   const std::vector<std::array<std::complex<double>, N>> &_data)
		: xgrid(_x), ygrid(_y), zgrid(_z), data(_data)
	{
	}

	/// Sets the grid and also sets data by a callable function oject
	InterpolationTable(std::vector<double> _x, std::vector<double> _y, std::vector<double> _z,
					   const std::function<std::array<std::complex<double>, N>(double x, double y,
																			   double z)> &function)
		: xgrid(_x), ygrid(_y), zgrid(_z)
	{
		// Set data using provided function
		data.resize(grid_size());
		int index = 0;
		for (int i = 0; i < xgrid.size(); i++)
		{
			double x = xgrid[i];
			for (int j = 0; j < ygrid.size(); j++)
			{
				double y = ygrid[j];
				for (int k = 0; k < zgrid.size(); k++)
				{
					double z = zgrid[k];
					data[index++] = function(x, y, z);
				}
			}
		}
	}

	/// Compute the function value by interpolation at the specified point
	std::array<std::complex<double>, N> compute_at(double x, double y, double z) const
	{
		Int3D stencil_ijk = locate_stencil(x, y, z);
		std::array<std::complex<double>, N> result;
		result.fill(0.0);
		int LX = stencil_size[0];
		int LY = stencil_size[1];
		int LZ = stencil_size[2];

		// This is a column vector, but the result should correspond to a row of the F matrix

		Int3D Lijk = {0, 0, 0};
		Double3D position = {x, y, z};
		for (int i = 0; i < LX; i++)
		{
			Lijk[0] = i;
			for (int j = 0; j < LY; j++)
			{
				Lijk[1] = j;
				for (int k = 0; k < LZ; k++)
				{
					Lijk[2] = k;
					Int3D data_ijk = Lijk + stencil_ijk;
					auto y_ijk = data[get_index(data_ijk)];
					double l_ijk = compute_lagrange_polynomial(stencil_ijk, Lijk, position);
					result += y_ijk * l_ijk;
				}
			}
		}
		return result;
	}

	/// Total number of grid points
	int grid_size() const { return xgrid.size() * ygrid.size() * zgrid.size(); }

	/// The number of stencils is directly related to the number of points.
	Int3D num_stencils() const
	{
		Int3D num_points = {xgrid.size(), ygrid.size(), zgrid.size()};
		return num_points - stencil_size + 1;
	}

	/// Helpful functions for converting from (i,j,k) indices to
	/// actual index in row major array, will check range in debug.
	int get_index(const Int3D &ijk) const
	{
		assert(ijk[0] >= 0 && ijk[1] >= 0 && ijk[2] >= 0);
		assert(ijk[0] < xgrid.size() && ijk[1] < ygrid.size() && ijk[2] < zgrid.size());
		const int ny = ygrid.size();
		const int nz = zgrid.size();
		return ijk[2] + nz * (ijk[1] + ny * ijk[0]);
	}

	/// Reverse of above, given an index returns the (i,j,k) indices, will check range in debug.
	Int3D get_ijk(const int index) const
	{
		assert(index >= 0);
		assert(index < grid_size());
		const int ny = ygrid.size();
		const int nz = zgrid.size();
		Int3D ijk;

		// Row Major - slightly optimized
		ijk[1] = index / nz;
		ijk[2] = index % nz;
		ijk[0] = ijk[1] / ny;
		ijk[1] = ijk[1] % ny;
		return ijk;
	}
	/// Number of interpolation points along each dimension. Equivalent to p+1, where p is the
	/// order of polynomial used for interpolation.
	Int3D stencil_size = {3, 3, 3};

  private:
	Int3D locate_stencil(double x, double y, double z) const
	{
		Int3D stencil_ijk;
		stencil_ijk[0] = locate_stencil_1d(x, xgrid, stencil_size[0]);
		stencil_ijk[1] = locate_stencil_1d(y, ygrid, stencil_size[1]);
		stencil_ijk[2] = locate_stencil_1d(z, zgrid, stencil_size[2]);
		return stencil_ijk;
	}

	static int locate_stencil_1d(double coord, std::vector<double> grid, int stencil_size)
	{
		auto begin = grid.cbegin();
		auto end = grid.cend();
		int stencil_1d = -1;
		int max_stencil = grid.size() - stencil_size;
		// Even stencils are centered on intervals
		if (stencil_size % 2 == 0)
		{
			stencil_1d = find_interval(coord, begin, end);
			// stencil indices are defined by the first point that belongs to them
			// So the width of the stencil is subtracted from the interval location
			stencil_1d -= stencil_size / 2 - 1;
		}
		else // Odd stencils are centered on grid points
		{
			stencil_1d = find_closest_point(coord, begin, end);
			stencil_1d -= (stencil_size - 1) / 2;
		}
		// Take care of stencil indices that were computed past the boundary
		stencil_1d = std::max(0, stencil_1d);
		stencil_1d = std::min(max_stencil, stencil_1d);

		return stencil_1d;
	}
	/// Computes one of the lagrange polynomials, given by Lijk, at a position within a given
	/// stencil.
	double compute_lagrange_polynomial(Int3D stencil_ijk, Int3D Lijk, Double3D pos) const
	{
		double product = 1.0;
		int LX = stencil_size[0];
		int LY = stencil_size[1];
		int LZ = stencil_size[2];
		// Compute product along x dimension
		for (int i = 0; i < LX; ++i)
		{
			if (i != Lijk[0])
			{
				const double x_m = xgrid[i + stencil_ijk[0]];
				const double x_idx = xgrid[Lijk[0] + stencil_ijk[0]];

				product *= (pos[0] - x_m) / (x_idx - x_m);
			}
		}
		// Compute product along y dimension
		for (int j = 0; j < LY; ++j)
		{
			if (j != Lijk[1])
			{
				const double y_m = ygrid[j + stencil_ijk[1]];
				const double y_idx = ygrid[Lijk[1] + stencil_ijk[1]];

				product *= (pos[1] - y_m) / (y_idx - y_m);
			}
		}
		// Compute product along z dimension
		for (int k = 0; k < LZ; ++k)
		{
			if (k != Lijk[2])
			{
				const double z_m = zgrid[k + stencil_ijk[2]];
				const double z_idx = zgrid[Lijk[2] + stencil_ijk[2]];

				product *= (pos[2] - z_m) / (z_idx - z_m);
			}
		}
		return product;
	}
	/// The 3D structured grid.
	std::vector<double> xgrid{}, ygrid{}, zgrid{};
	/// Stores many components or fields.
	std::vector<std::array<std::complex<double>, N>> data;
};

#endif // STRATA_INTERPOLATION_TABLE_HPP
