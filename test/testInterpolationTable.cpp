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

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "MGF.hpp"
#include "interpolation_table.hpp"

// Test polynomial should be perfectly interpolated by polynomial of order (1,2,3)
std::array<std::complex<double>, 1> polynomial(double x, double y, double z)
{
	return {0.5 * x + 0.125 * x * y * y + 1.6 * z * z * z + 3.1456 * x * y * y * z * z};
}

std::array<std::complex<double>, 1> green_hgf(double x, double y, double z)
{
	double k = 2 * strata::PI;
	double r = sqrt(x * x + y * y + z * z);
	std::complex<double> jkr = {0.0, k * r};
	return {std::exp(-jkr) / r};
}

void try_interpolation_green()
{
	int N = 10;
	int N_fine = 50 * N;
	double mid_rtol = 1e-4;
	std::vector<double> xgrid, ygrid, zgrid;
	double min = 0.1;
	double max = 2.0;
	strata::linspace(min, max, N, xgrid);
	strata::linspace(min, max, N, ygrid);
	strata::linspace(min, max, N, zgrid);

	InterpolationTable<1> tbl(xgrid, ygrid, zgrid, green_hgf);
	// Set stencil sizes that should perfectly interpolate data points
	tbl.stencil_size = {4, 4, 4};

	std::vector<double> x;
	strata::linspace(min, max, N_fine, x);
	double y = min;
	double z = min;
	double max_error = std::numeric_limits<double>::min();
	for (int i = 0; i < N_fine; i++)
	{
		auto reference = green_hgf(x[i], y, z);
		auto result = tbl.compute_at(x[i], y, z);
		double error = std::abs(reference[0] - result[0]) / std::abs(reference[0]);
		max_error = std::max(max_error, error);
	}
	std::cout << "Max Error: " << max_error << std::endl;

	bool refine = true;
	InterpolationTable<1> new_table;
	while(refine)
	{
		std::cout << "=====Try new table=====" << std::endl;
		max_error = std::numeric_limits<double>::min();
		refine = check_interpolation_and_update_grid<1>(tbl, new_table, green_hgf, mid_rtol);
		for (int i = 0; i < N_fine; i++)
		{
			auto reference = green_hgf(x[i], y, z);
			auto result = new_table.compute_at(x[i], y, z);
			double error = std::abs(reference[0] - result[0]) / std::abs(reference[0]);
			max_error = std::max(max_error, error);
		}
		std::cout << "Max Error: " << max_error << std::endl;
		tbl = new_table;
	}
}

int main(int argc, char **argv)
{
	try_interpolation_green();
	std::cout << "===========================" << std::endl;
	std::cout << "TestInterpolationTable()" << std::endl;
	std::cout << "===========================" << std::endl;

	std::string tech_file, out_file;

	if (argc < 2)
	{
		// throw std::runtime_error("[ERROR] Layer file was not provided.");
	}
	else if (argc < 3)
	{
		// tech_file = argv[1];
		// out_file = "MGFdata_interp.txt";
	}
	else
	{
		// tech_file = argv[1];
		// out_file = argv[2];
	}

	std::vector<double> xgrid{1, 2, 4, 5, 6}, ygrid{1, 2, 4, 5, 6}, zgrid{1, 2, 4, 5, 6};

	InterpolationTable<1> tbl(xgrid, ygrid, zgrid, polynomial);
	// Set stencil sizes that should perfectly interpolate data points
	tbl.stencil_size = {2, 3, 4};
	// Test point
	double x = 5.23541;
	double y = 3.8;
	double z = 1.89;

	auto reference = polynomial(x, y, z);
	auto result = tbl.compute_at(x, y, z);
	double error = std::abs(reference[0] - result[0]) / std::abs(reference[0]);
	std::cout << "Error: " << error << std::endl;

	if (error <= 10 * std::numeric_limits<double>::epsilon())
		return 0;
	else
		return 1;
}
