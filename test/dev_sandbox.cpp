// Author: Shashwat Sharma

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

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "MGF.hpp"

void save_mgf_table(const std::string &filename, MGF &mgf)
{
	std::ofstream outfile(filename);
	// For post-processing, we'll store the frequency and positions along the z axis in the header
	outfile << "# All SI units (Hz, m)" << std::endl;
	outfile << "frequency: " << mgf.f << std::endl;
	int layerid = (int)mgf.lm.z_nodes.size() - 1;
	int total_z_nodes = 0;
	//Save z nodes by layer
	outfile << "# Z nodes stored as (LayerID, z_coordinate)" << std::endl;
	for(int l_obs = (int)mgf.lm.z_nodes.size() - 1; l_obs >= 0; l_obs--)
	{
		total_z_nodes += (int)mgf.lm.z_nodes[l_obs].size();
		for(auto& znodes : mgf.lm.z_nodes[l_obs])
			outfile << layerid << " " << znodes << " ";
		layerid--;
	}
	outfile << std::endl;
	outfile << "total_z_nodes: " << total_z_nodes << std::endl;
	//Save rho nodes
	outfile << "rho_nodes: ";
	for(auto& rho : mgf.lm.rho_nodes)
		outfile << rho << " ";
	outfile << std::endl;
	outfile << "total_rho_nodes: " << (int)mgf.lm.rho_nodes.size() << std::endl;

	outfile << "# Data Table Format (z_observe, z_source, rho, component) stored as complex number row-major."  << std::endl;
	outfile << "# Green's function components (Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi)"  << std::endl;

	std::array<std::complex<double>, 9> G_dyadic;
	std::complex<double> G_phi;
	for (int l_obs = (int)mgf.lm.z_nodes.size() - 1; l_obs >= 0; l_obs--)
		for (int z_obs = 0; z_obs < mgf.lm.z_nodes[l_obs].size(); z_obs++)
		{
			double z_obs_pos = mgf.lm.z_nodes[l_obs][z_obs];
			for (int l_src = (int)mgf.lm.z_nodes.size() - 1; l_src >= 0; l_src--)
				for (int z_src = 0; z_src < mgf.lm.z_nodes[l_src].size(); z_src++)
				{
					double z_src_pos = mgf.lm.z_nodes[l_src][z_src];
					mgf.SetLayers(mgf.lm.FindLayer(z_src_pos), mgf.lm.FindLayer(z_obs_pos));
					//mgf.SetLayers(l_src, l_obs);
					for (auto &rho : mgf.lm.rho_nodes)
					{
						mgf.ComputeMGF(rho, 0.0, z_obs_pos, z_src_pos, G_dyadic, G_phi);
						for (int i = 0; i < 9; i++)
							outfile << G_dyadic[i].real() << " " << G_dyadic[i].imag() << " ";
						outfile << G_phi.real() <<  " " << G_phi.imag() << std::endl;
					}
				}
		}
}

int main(int argc, char** argv)
{

	std::cout << "===========================" << std::endl;
	std::cout << "dev_sandbox()" << std::endl;
	std::cout << "===========================" << std::endl;

    std::string tech_file, out_file;
	
	if (argc < 2)
	{
		throw std::runtime_error("[ERROR] Layer file was not provided.");
	}
	else if (argc < 3)
	{
		tech_file = argv[1];
		out_file = "MGFdata";
	}
	else
	{
		tech_file = argv[1];
		out_file = argv[2];
	}


	// ====== Input settings ======

	// ------ Ling, Jin, 2000 ------

	// // Set the analysis frequency and wave number
	// double f = 30.0e9;

	// double x_src = 0.0, y_src = 0.0, z_src = 0.4e-3;
	// double y_obs = 0.0, z_obs = 1.4e-3;
	
	// int Nx = 500; // Number of points in the sweep
	// double x_obs_min = std::abs(1.6e-4*lambda0);
	// double x_obs_max = std::abs(1.6e1*lambda0);

	// ------ AMD tcoil 2016 ------

	// Set the analysis frequency and wave number
	double f = 40e9;

	double x_src = 0.0, y_src = 0.0, z_src = -20.72e-6;
	double y_obs = 0.0, z_obs = 0.0;
	
	int Nx = 200; // Number of points in the sweep
	double x_obs_min = std::abs(5.0e-6);
	double x_obs_max = std::abs(16.0e-3);

	// ====== Stackup ======

	// Create a layer management object and parse the layer file
	LayerManager lm;
	lm.ProcessTechFile(tech_file, 1.0e-6);

	double omega = 2.0*M_PI*f;

	// Some useful constants are provided via the Strata namespace
	double k0 = omega*std::sqrt(strata::eps0*strata::mu0);
	double lambda0 = 2.0*M_PI/k0;

	// Precompute frequency-dependent layer data
	lm.ProcessLayers(f);

	// Print layer data to the terminal for verification
	lm.PrintLayerData(NULL, true);

	int i = lm.FindLayer(z_src);
	int m = lm.FindLayer(z_obs);

	// ====== Set up the source and observation points ======
	//std::vector<double> znodes_sample = {-38.22e-6, -29.47e-6, -20.72e-6, 0.0, 15.63e-6};
	//std::vector<double> znodes_sample = {-38.22e-6, -24.7575e-6, -11.295e-6, 2.1675e-6, 15.63e-6};
	std::vector<double> znodes_sample = {0.5e-6};
	lm.InsertNodes_z(znodes_sample);
	// We can use the Matlab-like linspace or logspace functions to create linearly- or logarithmically-spaced vectors points, provided via the Strata namespace
	std::vector<double> x_vec;
	//strata::logspace(std::log10(x_obs_min), std::log10(x_obs_max), Nx, x_vec);
	strata::linspace(x_obs_min, x_obs_max, Nx, x_vec);
//	if (i == m)
//	{
//		std::vector<double> z_nodes = {z_obs, z_src};
//		lm.InsertNodes_z(z_nodes, i);
//	}
//	else
//	{
//		std::vector<double> z_nodes = {z_src};
//		lm.InsertNodes_z(z_nodes, i);
//		z_nodes = {z_obs};
//		lm.InsertNodes_z(z_nodes, m);
//	}

	int N_rho = std::max(30.0, 30.0*x_obs_max/lambda0);
	std::vector<double> rho_nodes;
	strata::linspace(std::sqrt(std::pow(x_obs_min, 2) + std::pow(y_obs, 2)), std::sqrt(std::pow(x_obs_max, 2) + std::pow(y_obs, 2)), N_rho, rho_nodes);
	rho_nodes.insert(rho_nodes.cbegin(), 0.0);
	rho_nodes = x_vec;
	lm.ClearNodes_rho();
	lm.InsertNodes_rho(rho_nodes);
	lm.PrintNodeData_rho(nullptr,true);
	lm.PrintNodeData_z(nullptr,true);

	// ====== Initialize the MGF class ======

	// In this example, we'll compute the MGF using straightforward numerical integration
	
	// This class stores all the settings we want to use
	MGF_settings s;
	
	// This class is the "engine" which will compute the MGF
	MGF mgf;

	// Tell the class that we want to use the numerical integration method
	s.method = MGF_INTERPOLATE;
	s.extract_quasistatic = true;
	s.extract_singularities = false;
	s.extract_homogeneous = true;

	s.tol_qse = -1.0e-1;
	s.switching_point = -1.0;

	s.DCIM_method = DCIM_TWO_LEVEL;
	s.tol_svd = 1.0e-4;
	s.tol_eig = 1.0e-16;
	s.max_num_images = -1;

	s.sampling_method = MGF_INTEGRATE;
	s.order = 2;
	s.load_table = false;
	s.export_table = true;
	s.filename = out_file + ".mgf";

	// Initialize the class with the given stackup and chosen settings
	mgf.Initialize(f, lm, s);

	// The MGF class needs to know which layers we're working with, in order to perform some precomputations which may save time later on.
	// The working source and observation layers can be found with the FindLayer() method of the layer manager.
	// In a realistic MoM setting, if we are looping through the points on the mesh of an object, and we know in which layer that object resides, we can just set the source and observation layers once for each pair of source and observation objects.
	//int i = 0;// lm.FindLayer(z_src);
	//int m = 0;//lm.FindLayer(z_obs);
	mgf.SetLayers(i, m); // Source first, observation second


	// ====== Compute the MGF for each observation point ======
	
	// Create an output file where the MGF will be exported for post-processing
	std::string real_file = out_file + "_real.txt";
	std::ofstream outfile(real_file);

	// For post-processing, we'll store the frequency and positions along the z axis in the header
	outfile << "Frequency: " << f << " Hz" << std::endl;
	outfile << "z_src: " << z_src << " m" << std::endl;
	outfile << "z_obs: " << z_obs << " m" << std::endl;

	// The lateral separation between source and observation points (rho) and all the MGF components is tabulated for each observation point
	outfile << "\nrho Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;

	std::string imag_file = out_file + "_imag.txt";
	std::ofstream outfile_im(imag_file);

	// For post-processing, we'll store the frequency and positions along the z axis in the header
	outfile_im << "Frequency: " << f << " Hz" << std::endl;
	outfile_im << "z_src: " << z_src << " m" << std::endl;
	outfile_im << "z_obs: " << z_obs << " m" << std::endl;

	// The lateral separation between source and observation points (rho) and all the MGF components is tabulated for each observation point
	outfile_im << "\nrho Gxx Gxy Gxz Gyx Gyy Gyz Gzx Gzy Gzz Gphi" << std::endl;

	for (int ii = 0; ii < Nx; ii++)
	{

		double x_obs = x_vec[ii];

		// In the x and y directions, the MGF only depends on the separation between source and observation points, rather than the actual coordinates
		double x_diff = x_obs - x_src;
		double y_diff = y_obs - y_src;

		// The 3x3 dyadic MGF has 9 components which will be stored in a C++ standard array, in row major order.
		// In addition, there is a scalar component which is just a complex number.
		std::array<std::complex<double>, 9> G_dyadic;
		std::complex<double> G_phi;

		// Compute the MGF
		mgf.ComputeMGF(x_diff, y_diff, z_obs, z_src, G_dyadic, G_phi);

		// The dyadic MGF components are now stored in G_dyadic, while the scalar MGF is stored in G_phi.


		// ====== Optional modifications to the output ======

		// With reference to Michalski, Zheng, TAP 1990, 38 (3), equations (50)--(53), the MGF we computed above ** does include ** the cos and sin pre-factors. However, in literature, the MGF is often reported without these pre-factors. Therefore, for the purpose of this example, and to make direct comparisons to data from literature, those prefactors are cancelled out below. This section of code would not be needed in an actual MoM-type application.

//		std::complex<double> zeta = std::atan2(y_diff, x_diff);
//		std::complex<double> cos_term = std::cos(zeta);
//		std::complex<double> sin_term = std::sin(zeta);
//
//		if (std::abs(cos_term) > 0.0)
//		{
//			G_dyadic[2] /= cos_term;
//			G_dyadic[6] /= cos_term;
//		}
//		if (std::abs(sin_term) > 0.0)
//		{
//			G_dyadic[5] /= sin_term;
//			G_dyadic[7] /= sin_term;
//		}

		
		// ====== Export data to the output text file ======

		// Lateral separation
		double rho = std::sqrt( std::pow(x_diff, 2) + std::pow(y_diff, 2) );
//		outfile << rho << " " <<
//				std::abs(G_dyadic[0]) << " " << std::abs(G_dyadic[1]) << " " << std::abs(G_dyadic[2]) << " " <<
//				std::abs(G_dyadic[3]) << " " << std::abs(G_dyadic[4]) << " " << std::abs(G_dyadic[5]) << " " <<
//				std::abs(G_dyadic[6]) << " " << std::abs(G_dyadic[7]) << " " << std::abs(G_dyadic[8]) << " " << std::abs(G_phi) << std::endl;
		outfile << rho << " " <<
			std::real(G_dyadic[0]) << " " << std::real(G_dyadic[1]) << " " << std::real(G_dyadic[2]) << " " <<
			std::real(G_dyadic[3]) << " " << std::real(G_dyadic[4]) << " " << std::real(G_dyadic[5]) << " " <<
			std::real(G_dyadic[6]) << " " << std::real(G_dyadic[7]) << " " << std::real(G_dyadic[8]) << " " << std::real(G_phi) << std::endl;

		outfile_im << rho << " " <<
			std::imag(G_dyadic[0]) << " " << std::imag(G_dyadic[1]) << " " << std::imag(G_dyadic[2]) << " " <<
			std::imag(G_dyadic[3]) << " " << std::imag(G_dyadic[4]) << " " << std::imag(G_dyadic[5]) << " " <<
			std::imag(G_dyadic[6]) << " " << std::imag(G_dyadic[7]) << " " << std::imag(G_dyadic[8]) << " " << std::imag(G_phi) << std::endl;
	}

	save_mgf_table(out_file + "_complex.txt", mgf);

	return 0;
}






