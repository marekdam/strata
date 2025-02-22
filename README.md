Strata is a C++ library for computing the dyadic multilayer Green's function (MGF), intended for use in computational electromagnetics codes for solving Maxwell's equations in integral form.
Strata can be used as a standalone library to compute the MGF as a function of frequency and/or spatial separation.
It can also be incorporated easily into existing integral equation solvers.

## Prerequisites

* [CMake](https://cmake.org/), at least version 3.10.2

## Installation

To install Strata, clone this repository and navigate to the root directory of the cloned repository in a terminal. Then run the following commands:

```bash
mkdir build
cd build
cmake ..
make
```

**Important notes**

* An internet connection is required during the installation process, because some third-party libraries ([OpenBLAS](https://www.openblas.net/) and [yaml-cpp](https://github.com/jbeder/yaml-cpp)) will automatically be downloaded and installated locally in the directory `external`.
* If you already have OpenBLAS installed on your machine, you could save some time by reusing it. In that case, instead of `cmake ..`, run
   ```bash
   cmake -DWITH_OPENBLAS_INC=location/of/openblas/headers -DWITH_OPENBLAS_LIB=location/of/openblas/lib ..
   ```
   where `location/of/openblas/headers` is the path to the directory containing all OpenBLAS-related headers (e.g., `cblas.h` and `lapacke.h`), and `location/of/openblas/lib` is the path to the directory containing the OpenBLAS shared library, usually `libopenblas.so`.
* On macOS, the `cmake` command requires an additional option:
   ```bash
   cmake -DCMAKE_INSTALL_RPATH_USE_LINK_PATH="ON" ..
   ```

Advanced users may want to consult [this page](https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#what-is-an-out-of-source-build) for out-of-source builds.

To use Strata within your own application, include the directory `inc`, which contains all the headers, and link to `build/libstrata.so` while building your project.

Strata is designed for UNIX-based operating systems, and has been tested on
* Ubuntu 20.04
* Ubuntu 18.04
* CentOS 7
* macOS Catalina 10.15.6

## Usage

Please consult the Strata documentation, `doc/strata.html` or `doc/strata.pdf` for usage instructions.

**Note**: the documentation is generated automatically during the build process.

### Run your first test case

A set of example layer definition files and test cases is provided in the `test` directory.
The most basic test case is in `test/testMGF.cpp`, which takes as input the layer definition file (see `doc/strata.html` or `doc/strata.pdf` for details) and a name for the output file.
From within the `build` folder, run

```bash
./testMGF ../test/examples/ling_jin_2000/layers.yaml ../test/examples/ling_jin_2000/MGFdata.txt
```

This will compute the MGF for a simple example and store the reuslts in the file `test/examples/ling_jin_2000/MGFdata.txt`.
Now, assuming you have Python 3 installed, you can plot the computed MGF:

```bash
cd ../test/examples/ling_jin_2000
python3 makeplots.py
```
You should see the following plot:

![Results of testMGF.cpp](doc/source/figures/plots.png?raw=true)

## Features

The functionality of Strata can be divided into the computation of the MGF in (a) the spectral domain, and (b) the spatial domain.

### a) Computation of the spectral-domain MGF

   Strata can compute, in spectral domain:

   * Each component of the dyadic **G<sup>(A)</sup>**, and the scalar term **G<sup>(&#966;)</sup>**, as defined in formulation-C of [[1]](#mgf01)
   * Each unique component of the dyadics **G<sup>(EM)</sup>** and **G<sup>(HJ)</sup>**, as defined in [[2]](#mgf02)

   In addition, one can optionally extract the quasistatic terms as defined in [[3]](#qse01).

### b) Computation of the spatial-domain MGF

   Strata allows computing the spatial domain MGF via the following methods:

   * Direct numerical integration of the Sommerfeld integrals, as described in [[4]](#mgf03)
   * The discrete complex image method (DCIM) [[5]](#dcim01)
   * Precomputation of interpolation tables constructed using either one of the above methods

   In addition, one can optionally extract (in spectral domain) and add (in spatial domain) the quasistatic terms as defined in [[3]](#qse01).

   To easily incorporate into existing codes, one can also extract out the singular behaviour of the MGF, so that the singularities can be treated separately.

Strata is written with modularity and readability as a priority, to allow users to incorporate new formulations and extend the functionality in both spectral and spatial domains.

## Citing Strata

We request that you acknowledge the authors of Strata by citing the following:

```latex
@INPROCEEDINGS{strata,
	author={S. {Sharma} and P. {Triverio}},
	booktitle={2021 {IEEE} International Symposium on Antennas and Propagation and {USNC-URSI} Radio Science Meeting},
	title={Strata: An Open-Source {C++} Library for Computing {Green's} Functions for Layered Media},
	year={2021},
	month={Dec.},
	address = {Singapore}}
```
## Issues and Improvements

We encourage users to report problems and suggest improvements by opening issues within the GitHub environment.

We also encourage users to help in our goal of providing a high-quality library for layered medium Green's functions by contributing features and new formulations.

## Known Pitfalls

* This library was originally developed for chip-level structures where conductors are embedded in dielectric layers, as opposed to printed circuit-like structures, where the metal is placed at the top of a substrate. I.e., there cannot be source or observation points in the upper half-space - all points must lie inside a layer. To model printed circuit structures, create a phantom air layer above the substrate, and place the structure within it. Note that no two layers are allowed to have identical material properties - typically, if this happens, the layers are automatically merged into a single layer. However, if both the phantom layer and the upper half space are supposed to be air, you must make the relative permittivity of the phantom layer slightly different from that of air (e.g., 1.001) to avoid numerical issues.

* Again, because the upper half-space is not considered a layer, our indexing scheme is offset from the scheme in [[1]](#mgf01). In [[1]](#mgf01), the indexing starts from the upper half-space, which has index 0. In Strata, the first layer has index 0, not the upper half-space.

* These comments also apply to the lower half-space.

## Related Projects

This is a list of related open-source projects by other authors - please let us know if some projects were missed.

* A Matlab and Fortran library for computing the MGF for the case of two half spaces can be found [here](https://github.com/UniPD-DII-ETCOMP/Half_Space_Green_A_Phi).

## References

<a name="mgf01"></a>[1] K. A. Michalski and D. Zheng, "Electromagnetic scattering and radiation by surfaces of arbitrary shape in layered media. I. theory," *IEEE Trans. Antennas Propag.*, vol. 38, no. 3, pp. 335–344, Mar. 1990.

<a name="mgf02"></a>[2] K. A. Michalski and J. R. Mosig, "Multilayered media Green's functions in integral equation formulations," *IEEE Trans. Antennas Propag.*, vol. 45, no. 3, pp. 508–519, Mar. 1997.

<a name="qse01"></a>[3] E. Simsek, Q. H. Liu, and B. Wei, "Singularity Subtraction for Evaluation of Green's Functions for Multilayer Media," *IEEE Trans. Microw. Theory Tech.*, vol. 54, pp. 216–225, Jan 2006.

<a name="mgf03"></a>[4] K. A. Michalski and J. R. Mosig, "Efficient computation of Sommerfeld integral tails - methods and algorithms," *J. Electromagn. Waves Appl.*, vol. 30, no. 3, pp. 281–317, 2016.

<a name="dcim01"></a>[5] M. I. Aksun, "A Robust Approach for the derivation of closed-form Green's functions," *IEEE Trans. Microw. Theory Tech.*, vol. 44, pp. 651–658, May 1996.

## Contributors

Strata was developed as part of a PhD project at the [Modelics Lab](http://modelics.org/), University of Toronto, Canada.

* Shashwat Sharma (Architect & Developer)
* Piero Triverio (Principal Investigator)

