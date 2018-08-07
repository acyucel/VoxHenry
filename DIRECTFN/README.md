# DIRECTFN 
---

DIRECTFN package is free software for the accurate and efficient evaluation of 4-D singular integrals arising from Galerkin MoM surface integral 
equation formulations over conforming triangular or quadrilateral meshes.  The fully-numerical algorithms of DIRECTFN are suitable for the following applications: 

* Weakly and strongly singular kernels
* Planar and curvilinear elements
* Basis/testing functions of arbitrary order
* Problem-specific Green's functions (e.g. expressed in spectral integral form)
* From zero frequency to frequencies beyond microwaves
* Spectral convergence to machine precision

The complete manual for the use of the library can be found in [`ReadMe.pdf`](/docs/Manual).
## Download 

The source code of DIRECTFN is written in C++11. Also, mex plugins are available to provide a fast Matlab interface to this functionality. 
You are welcome to report any bugs found while testing!

This repository is organized as follows:

* [`./include`](/include)
folder contains headers of the library with declarations only.
*  [`./src`](/src) folder contains source files where all the templates are instantiated in corresponding files properly. So, in the client programs, one can use the templated classes which have been precompiled already. 
* [`./lib`](/lib) is the target path to store the compiled static library, which can be further used, e.g. in Matlab interface, or Python wrappers.
* [`./tests`](/tests) directory contains some
initial examples how to use the library in your own code. 
* [`./examples`](/examples) folder contains the scripts to reproduce the results, presented in paper [2].

## Compilation

The provided static library and examples can be compiled using [make](https://www.gnu.org/software/make/) utility. Therefore, it requires Windows users to have some Linux-like environment, such as [Cygwin](https://www.cygwin.com/), or [MinGW](http://www.mingw.org/) and [MSYS2](https://www.msys2.org/). However, the other possible option is to use [nmake](https://msdn.microsoft.com/en-us/library/dd9y37ha.aspx) utility from Visual Studio (see [`Compilation with nmake`](#compilation-with-nmake)). 

The following steps are necessary for compilation:

* If you use MSYS2, launch the MSYS2 console, and add to the path your compiler environment. If you are using Octave, you must use the same compiler that you will use to generate the Mex files in Octave. In this case, add to the path the Octave-provided version of GCC, e.g. "<OctaveInstallationPath>\Octave-4.4.0\bin" where <OctaveInstallationPath> is your Octave installation path, and assuming Octave version 4.4.0.

* Go to the [`./settings`](/settings) directory and create or copy-modify platform-specific `Makefile.your_cpu_alias`. Here the `your_cpu_alias` is an alias for your system compiler (like intel, gcc620, clang or so). This alias is just an abbreviation to distinguish different Makefiles and can be any name you like.

* Define environment variable `CPU=your_cpu_alias`. It can be done by the command
```
 export CPU=your_cpu_alias
```
or by figuring out the CPU explicitly when compiling the code with make utility:
```
 make CPU=your_cpu_alias
```

* Next, add  the following lines to the end of Makefile.in:

```
ifeq ($(CPU),your_cpu_alias)
    include $(MF_PREFIX)/Makefile.target
endif
```

* Finally, go to the [`./lib/unix`](github/lib/unix) folder and type 
```
 make
``` 
or 
```
 make CPU=your_cpu_alias
```
This command will compile the code and move the static library into the [`./lib`](/lib) folder.
Default name for the library is `directfn` (thus, it will be compiled and 
built in the static library file `libdirectfn.a`).
You can change the name in `./settings/Makefile.in` by defining
`DIRECTFN_LIB_NAME` variable.

# compilation-with-nmake

In order to enable working with the [nmake](https://msdn.microsoft.com/en-us/library/dd9y37ha.aspx) utility, we provide also special makefiles for it, since the syntax is slightly different from the one that [make](https://www.gnu.org/software/make/) uses. To compile the source code with nmake you should do the following simple steps:

 * The folder [`./lib/win/`](/lib/win) contains the Visual Studio project file to build the static library.
 * To compile the examples run the [Developer Command Prompt for Visual Studio](https://msdn.microsoft.com/en-us/library/ms229859(v=vs.110).aspx) or Visual Studio x86(x64) Native Tools, change the folder to `.examples/<Paper Name>/c++` and type

```
nmake -f Makefile_win
```

Once the library is built, you may go to the [`./tests`](/tests) or [`./examples`](/examples) folders to compile and run some examples. 

## Library Usage

The complete manual of how to use the existing code is presented in [`./docs/Manual.pdf`](/docs/Manual/). The basic patterns of using DIRECTFN library can be found in the folder [`./tests/basic`](github/tests/basic). Here we briefly describe the main algorithm of computing the surface-surface singular integrals. 
In C++, you first `#include "directfn_quad.h"`,  if you need to work with quadrilateral elements, or `#include "directfn_triag.h"`, if you need to work with triangles. Then you should define a "singular contour" `SingularContour3xn` according to the specific adjacency type with the proper number of points:
```
SingularContour3xn cntr;
cntr.set_points(r1, r2, r3, r4);
```
Here `double[3] ri` are coordinates of vertices of the considered elements in 3D--space.
 Note that the order of vertices for each case of adjacency must correspond to the rules defined below.
 Next, create an object of ST, EA or VA algorithm and use its interface instantiated by a `ParticularKernel` kernel type. 
```
unique_ptr<Quadrilateral_ST<ParticularKernel>> up_quad_st(new Quadrilateral_ST<ParticularKernel>());
```
The type `ParticularKernel`
is a C++ type which inherits the `AbstractKernel` class defined in `directfn_kernel_base.h/cpp`. The specific kernels implemented for quadrilaterals and triangles are discussed in details in the next sections.
Now you should define the wavenumber `k0` and the orders of Gauss-Legendre quadratures `N1, N2, N3, N4`, and assign your contour `cntr` with vertices to your object `up_quad_st`.
```
up_quad_st.set_wavenumber(k0);
up_quad_st.set_Gaussian_orders_4(N1, N2, N3, N4);
up_quad_st.set(cntr);
```
Now you are ready to calculate the surface-surface singular integral(s) by calling
```
up_quad_st.calc_Iss(); // Iss = Integral Surface-Surface
```
and use the obtained values.
```
const dcomplex * ref_val = up_quad_st->Iss();
// use ref_val somewhere
```

In case you need to recalculate these integrals for several contours or wave numbers we
strongly recommend to do it as it is shown above: create an object of the
`Quadrilateral_ST(EA,VA)` algorithm once and then use it for contour points or other parameters
reset.  This helps to avoid time-consuming memory allocation/deallocation related to creating and destroying the object.
In case you need to compute the integral once,
you can use a simpler interface, developed for Matlab scripts.

## Examples

The folder [`./examples`](/examples) contains the scripts to reproduce results, presented in paper [2]. For a given example there are two ways to execute it:

* Go to [`./examples/EuCap2017_Examples/c++`](/examples/EuCap2017_Examples/c++) or [`./examples/Full_Paper_Examples/c++] (/examples/Full_Paper_Examples/c++) folder and compile the C++ code:
```
make
```
or
```
nmake -f Makefile_win
```
if you are using nmake utility.
then type the name of the executable file
```
./test_name
``` 
It will create a file `Results_<test_name>.txt` with timings and errors, and you can  open `Examples.ipynb` notebook and visualize the results by executing the corresponding cell. Note that to run the .ipynb file you need [ipython](https://ipython.org/) to be installed.

* The other way is to use Matlab scripts. To do this, you first need to build the mex interface to C++ functions by going to the [`./mex`](/mex) folder and executing the `build.m` script. After that, you can go to the folder [`./examples/EuCap2017_Examples/matlab/`](/examples/EuCap2017_Examples/matlab/) or [`./examples/Full_Paper_Examples/matlab/`](/examples/Full_Paper_Examples/matlab/) and run any `test_name.m` script and it will do all the computations and produce the error plots and timings automatically.
* Note that if you're going to build mex files with one of the standart Windows compiler, e.g., provided by Visual Studio, you should build the static library with Windows compiler as well, i.e., using [`nmake`](https://msdn.microsoft.com/en-us/library/dd9y37ha.aspx).



## References

[1] A. G. Polimeridis, F. Vipiana, J. R. Mosig, and D. R. Wilton, “DIRECTFN:
Fully numerical algorithms for high precision computation
of singular integrals in Galerkin SIE methods,” *IEEE Trans. Antennas
Propag.*, vol. 61, no. 6, pp. 3112–3122, Jun. 2013.

[2] A. A. Tambova, M. S. Litsarev, G. D. Guryev, and A. G. Polimeridis, "On the generalization of DIRECTFN for singular integrals over quadrilateral patches", *IEEE Trans. Antennas
Propag.*, (submitted), [arXiv:1703.08146 [physics.comp-ph]](https://arxiv.org/abs/1703.08146).

## License

Copyright © 2016 Athanasios Polimeridis

DIRECTFN is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License (LGPL) as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

DIRECTFN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU GPLv3 for more details.

You should have received a copy of the GNU Lesser General Public License along with this program. If not, see http://www.gnu.org/licenses/.
