==================================================================================
   VoxHenry is an Inductance Extraction Simulator for Voxelized Geometries
           by Abdulkadir C. Yucel and Jacob K. White (MIT).
        Empowered by modules/ideas from Athanasios G. Polimeridis,
          Ioannis P. Georgakis (Skoltech), Hakan Bagci (KAUST),
and optimized by E. Di Lorenzo (FastFieldSolvers S.R.L., www.fastfieldsolvers.com)
==================================================================================

---------------------------
How-to-compile instructions
---------------------------

VoxHenry is written for MatLab, and has been tested under version 2018a,
both under Windows and under Linux.
VoxHenry has also been ported under Octave 4.4.0, but this version is slightly
less optimized, as it does not use LTDT decomposition for the preconditioner
and cannot use the parallelizaiton enabled by the 'parfor' in Matlab,
if the MatLab 'parallel toolbox' is available.

VoxHenry uses two additional modules (included in the distribution):

- DIRECTFN
- SuiteSparse

DIRECTFN
--------

DIRECTFN is not needed any more in the default VoxHenry version,
as the values needed from DIRECTFN are pre-calculated.

If you want to use DIRECTFN on the fly, be aware that DIRECTFN
comes in two flavors:

- one is faster, and aligned to the latest DIRECTFN distribution available
  from https://github.com/thanospol/DIRECTFN as of Aug 2018, but it requires
  compilation in format of a library before being able to use it
  from MatLab or Octave.
  
- one is slower, but can be directly compiled from within MatLab or Octave.
  If you want to use this version, you must modify the two files:
  - compile_mex_routines.m
  - pre_define_the_path_for_folders.m
  changing "use_recent_DIRECTFN=1" in "use_recent_DIRECTFN=0"

If you are using the fast DIRECTFN implementation, you need to compile it
first in the form of a library. Follow the instructions you can find in the
file "README.md" under the directory "DIRECTFN" to compile the library
for your environment.

Then, before using VoxHenry, you should launch the script 'compile_mex_routines.m'
contained in the root directory, to generate the corresponding 'mex' files
(i.e. C/C++ compiled files that can be called from MatLab/Octave).

SuiteSparse
-----------

SuiteSparse is distributed in compiled format for Windows and Linux,
usable by MatLab. SuiteSparse is currenty not used in the Octave port.


You are welcome to download and compile more recent versions of the packages.

Please check the license conditions of these modules if you are interested
in redistribution.

---------------------------
Usage instructions
---------------------------

To launch a simulation, run the file

VoxHenry_executer.m

and select the input file you want to use for the simulation from the list
that will be displayed on the screen.
The input files listed are contained in the sub-directory 'Input_files',
with the extension '.vhr'
You can add your input files to this sub-directory. The format of the .vhr
files is pretty simple and is described by the comments contained in any of the 
default input files. You can manually create your own input file, 
or you can use any of the four generator scripts contained in the VoxHenry
distribution to create variations of the input files.
The generators are:

- VoxHenry_generator_numex1_straight_conductor.m
- VoxHenry_generator_numex2_wire.m
- VoxHenry_generator_numex3_circular_coil.m
- VoxHenry_generator_numex4_square_coil.m

You can change the parameters in the relevant section of these files
(e.g. dimension of the voxels, frequency range, number or conductors, etc.)
or use them as the basis for your own geometry creation.

The results are printed on the screen, as well as stored in the format
of MatLab matrix files under the sub-directory 'Results'.

Note that by default VoxHenry_executer.m plots a drawing of the
current densities at one defined frequency. This operation may be
time consuming for large simulations. You can avoid that by changing
the relevant switch in the VoxHenry_executer.m input file.

