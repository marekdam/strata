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
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "MGF.hpp"
#include "interpolation_table.hpp"

void save_to_file(std::vector<double> &x, std::vector<double> &y, std::string filename)
{
	std::ofstream file(filename);
	for (int i = 0; i < x.size(); i++)
	{
		file << x[i] << "\t";
	}
	file << std::endl;
	for (int i = 0; i < y.size(); i++)
	{
		file << y[i] << "\t";
	}
	file.close();
}

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

int try_interpolation_green()
{
	int N = 10;
	double mid_rtol = 1.0e-4;
	std::vector<double> xgrid, ygrid, zgrid;
	double min = 0.0001;
	double max = 2;
	strata::linspace(min, max, N, xgrid);
	strata::linspace(min, max, N, ygrid);
	strata::linspace(min, max, N, zgrid);

	InterpolationTable<1> tbl(xgrid, ygrid, zgrid, green_hgf);
	// Set stencil sizes that should perfectly interpolate data points
	tbl.stencil_size = {4, 4, 4};

	double y = min;
	double z = min;

	bool refine = true;
	InterpolationTable<1> new_table;
	new_table.stencil_size = {4, 4, 4};
	int tbl_iter = 1;
	while (refine)
	{
		std::cout << "Table: " << tbl_iter++;
		double max_error = 0.0;
		refine =
			check_interpolation_and_update_grid<1>(tbl, new_table, green_hgf, max_error, mid_rtol);
		std::cout << "\tMax error encountered: " << max_error << std::endl;
		tbl = new_table;
	}

	std::vector<double> xvec, yvec, zvec;
	tbl.get_grids(xvec, yvec, zvec);
	std::vector<double> ref_result(2 * (xvec.size() - 1)), tbl_result(2 * (xvec.size() - 1));
	std::vector<double> ref_arg_result(2 * (xvec.size() - 1)),
		tbl_arg_result(2 * (xvec.size() - 1));
	std::vector<double> x_midpoints(2 * (xvec.size() - 1));
	for (int i = 0; i < xvec.size() - 1; i++)
	{
		double xdelta = 0.5 * (xvec[i + 1] - xvec[i]);
		x_midpoints[2 * i] = xvec[i];
		ref_result[2 * i] = std::abs(green_hgf(xvec[i], y, z)[0]);
		tbl_result[2 * i] = std::abs(new_table.compute_at(xvec[i], y, z)[0]);
		ref_arg_result[2 * i] = std::arg(green_hgf(xvec[i], y, z)[0]);
		tbl_arg_result[2 * i] = std::arg(new_table.compute_at(xvec[i], y, z)[0]);
		x_midpoints[2 * i + 1] = xvec[i] + xdelta;
		ref_result[2 * i + 1] = std::abs(green_hgf(xvec[i] + xdelta, y, z)[0]);
		tbl_result[2 * i + 1] = std::abs(new_table.compute_at(xvec[i] + xdelta, y, z)[0]);
		ref_arg_result[2 * i + 1] = std::arg(green_hgf(xvec[i] + xdelta, y, z)[0]);
		tbl_arg_result[2 * i + 1] = std::arg(new_table.compute_at(xvec[i] + xdelta, y, z)[0]);
	}
	save_to_file(x_midpoints, ref_result, "ref.txt");
	save_to_file(x_midpoints, tbl_result, "tbl.txt");
	save_to_file(x_midpoints, ref_arg_result, "ref_arg.txt");
	save_to_file(x_midpoints, tbl_arg_result, "tbl_arg.txt");

	tbl.get_grids(xvec, yvec, zvec);
	std::vector<double> count(xvec.size() - 1), xdelta(xvec.size() - 1);
	std::iota(count.begin(), count.end(), 0);
	for (int i = 0; i < xvec.size() - 1; i++)
	{
		xdelta[i] = xvec[i + 1] - xvec[i];
	}

	save_to_file(count, xdelta, "interpolation.txt");
	return 0;
}

int test_polynomial_interpolation()
{
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

int test_MGF_interpolation()
{
	std::string tech_file = TEST_DIR;
	tech_file += "ling_jin_2000/layers.yaml";
	LayerManager lm;
	lm.ProcessTechFile(tech_file);

	// Set the analysis frequency and wave number
	double f = 30.0e9;
	double omega = 2.0 * M_PI * f;

	// Some useful constants are provided via the Strata namespace
	double k0 = omega * std::sqrt(strata::eps0 * strata::mu0);
	double lambda0 = 2.0 * M_PI / k0;

	lm.ProcessLayers(f);
	lm.PrintLayerData(NULL, true);

	// ------ Ling, Jin, 2000 ------

	// For this example, we'll sweep the observation point along the x axis from 10^{-4} wavelengths
	// to 10 wavelengths away from the source point

	double x_src = 0.0, y_src = 0.0, z_src = 0.4e-3;
	double y_obs = 0.0, z_obs = 1.4e-3;
	z_obs = 1.1*z_src;
	int Nx = 500; // Number of points in the sweep
	//double x_obs_min = std::abs(1.6e-4 * lambda0);
	double x_obs_min = std::abs(1.6e-4 * lambda0);
	//double x_obs_max = std::abs(1.6e1 * lambda0);
	double x_obs_max = std::abs(1.0e0 * lambda0);

	// ------ Generate the points ------

	// We can use the Matlab-like linspace or logspace functions to create linearly- or
	// logarithmically-spaced vectors points, provided via the Strata namespace
	std::vector<double> x_vec;
	strata::logspace(std::log10(x_obs_min), std::log10(x_obs_max), Nx, x_vec);

	MGF_settings s;
	MGF mgf;
	s.method = MGF_INTEGRATE;
	s.extract_quasistatic = true;
	s.extract_homogeneous = true;
	s.tol_qse = 1e-6;
	mgf.Initialize(f, lm, s);

	int i = lm.FindLayer(z_src);
	int m = lm.FindLayer(z_obs);
	mgf.SetLayers(i, m); // Source first, observation second

	//double x_obs = x_vec[13];

	// In the x and y directions, the MGF only depends on the separation between source and observation points, rather than the actual coordinates
	//double x_diff = x_obs - x_src;
	//double y_diff = y_obs - y_src;

	//std::array<std::complex<double>, 9> G_dyadic;
	//std::complex<double> G_phi;

	// Compute the MGF
	//mgf.ComputeMGF(x_diff, y_diff, z_obs, z_src, G_dyadic, G_phi);

	auto mgf_wrap = [&mgf](double z_obs, double z_src, double rho)->std::array<std::complex<double>, 10>
	{
		std::array<std::complex<double>, 10> result;
		std::array<std::complex<double>, 9> G_dyadic;
		std::complex<double> G_phi;
		mgf.ComputeMGF(rho, 0, z_obs, z_src, G_dyadic, G_phi);
		std::copy_n(G_dyadic.cbegin(),9,result.begin());
		result[9] = G_phi;
		return result;
	};
	std::function<std::array<std::complex<double>, 10>(double, double, double)> mgf_func_wrap(mgf_wrap);

	int start_interpolation_N = 20;
	std::vector<double> zo_grid, zp_grid, rho_grid;
	strata::linspace(z_obs, 1.001*z_obs, 2, zo_grid);
	strata::linspace(z_src, 1.001*z_src, 2, zp_grid);
	strata::linspace(x_obs_min, x_obs_max, start_interpolation_N, rho_grid);

	InterpolationTable<10> tbl(zo_grid, zp_grid, rho_grid, mgf_func_wrap);
	tbl.stencil_size = {2, 2, 4};

	int tbl_iter = 1;
	bool refine = true;
	double mid_rtol = 0.001;
	InterpolationTable<10> new_table;
	new_table.stencil_size = tbl.stencil_size;
	while (refine && tbl_iter <= 10)
	{
		std::cout << "Table: " << tbl_iter++;
		double max_error = 0.0;
		refine =
			check_interpolation_and_update_grid<10>(tbl, new_table, mgf_func_wrap, max_error, mid_rtol);
		std::cout << "\tMax error encountered: " << max_error << std::endl;
		tbl = new_table;
	}

	std::vector<double> xvec, yvec, zvec;
	tbl.get_grids(xvec, yvec, zvec);
	std::vector<double> ref_result(2 * (zvec.size() - 1)), tbl_result(2 * (zvec.size() - 1));
	std::vector<double> ref_arg_result(2 * (zvec.size() - 1)),
		tbl_arg_result(2 * (zvec.size() - 1));
	std::vector<double> z_midpoints(2 * (zvec.size() - 1));
	int comp = 9;
	for (int i = 0; i < zvec.size() - 1; i++)
	{
		double zdelta = 0.5 * (zvec[i + 1] - zvec[i]);
		z_midpoints[2 * i] = zvec[i];
		ref_result[2 * i] = std::real(mgf_func_wrap(z_obs, z_src, zvec[i])[comp]);
		tbl_result[2 * i] = std::real(new_table.compute_at(z_obs, z_src, zvec[i])[comp]);
		ref_arg_result[2 * i] = std::imag(mgf_func_wrap(z_obs, z_src, zvec[i])[comp]);
		tbl_arg_result[2 * i] = std::imag(new_table.compute_at(z_obs, z_src, zvec[i])[comp]);
		z_midpoints[2 * i + 1] = zvec[i] + zdelta;
		ref_result[2 * i + 1] = std::real(mgf_func_wrap(z_obs, z_src, zvec[i] + zdelta)[comp]);
		tbl_result[2 * i + 1] = std::real(new_table.compute_at(z_obs, z_src, zvec[i] + zdelta)[comp]);
		ref_arg_result[2 * i + 1] = std::imag(mgf_func_wrap(z_obs, z_src, zvec[i] + zdelta)[comp]);
		tbl_arg_result[2 * i + 1] = std::imag(new_table.compute_at(z_obs, z_src, zvec[i] + zdelta)[comp]);
	}
	save_to_file(z_midpoints, ref_result, "ref.txt");
	save_to_file(z_midpoints, tbl_result, "tbl.txt");
	save_to_file(z_midpoints, ref_arg_result, "ref_arg.txt");
	save_to_file(z_midpoints, tbl_arg_result, "tbl_arg.txt");

	tbl.get_grids(xvec, yvec, zvec);
	std::vector<double> count(zvec.size() - 1), zdelta(zvec.size() - 1);
	std::iota(count.begin(), count.end(), 0);
	for (int i = 0; i < zvec.size() - 1; i++)
	{
		zdelta[i] = zvec[i + 1] - zvec[i];
	}

	save_to_file(count, zdelta, "interpolation.txt");

	return 0;
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

	return test_MGF_interpolation();
}
