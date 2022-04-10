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

//Test polynomial should be perfectly interpolated by polynomial of order (1,2,3)
std::array<std::complex<double>, 1> polynomial(double x, double y, double z)
{
	return {0.5 * x + 0.125*x * y * y + 1.6*z * z * z + 3.1456*x * y * y * z * z};
}

int main(int argc, char **argv)
{

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

	if (error <= 10*std::numeric_limits<double>::epsilon())
		return 0;
	else
		return 1;
}
