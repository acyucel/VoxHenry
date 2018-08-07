==================================================================================
VoxHenry is an Inductance Extraction Simulator for Voxelized Geometries
by Abdulkadir C. Yucel and Jacob K. White (MIT).
Empowered by modules/ideas from Athanasios G. Polimeridis and Ioannis P. Georgakis (Skoltech), 
Hakan Bagci (KAUST), and optimized by E. Di Lorenzo (www.fastfieldsolvers.com)
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

However, DIRECTFN comes in two flavors:

- one is faster, and aligned to the latest DIRECTFN distribution available
  from https://github.com/thanospol/DIRECTFN as of Aug 2018, but it requires
  compilation in format of a library before being able to use it
  from MatLab or Octave. This is the default one.
  
- one is slower, but can be directly compiled from within MatLab or Octave.
  If you want to use this version, you must modify the two files:
  - compile_mex_routines.m
  - pre_define_the_path_for_folders.m
  changing "use_recent_DIRECTFN=1" in "use_recent_DIRECTFN=0"

If you are using the default DIRECTFN implementation, you need to compile it
first in the form of a library. Follow the instructions you can find in the
file "README.md" under the directory "DIRECTFN" to compile the library
for your environment.

Then, before using VoxHenry, you should launch the script 'compile_mex_routines.m'
contained in the root directory, to generate the corresponding 'mex' files
(i.e. C/C++ compiled files that can be called from MatLab/Octave).

You are welcome to download and compile more recent versions of the packages.

Please check the license conditions of these modules if you are interested
in redistribution.

---------------------------
Usage instructions
---------------------------

You can run any of the four main files 

- VoxHenry_executer_numex1_straight_conductor.m
- VoxHenry_executer_numex2_wire.m
- VoxHenry_executer_numex3_circular_coil.m
- VoxHenry_executer_numex4_square_coil.m

to analyze the corresponding geometries.
Open these files to enable / disable specific simulation options
(e.g. dimension of the voxels, frequency range, number or conductors, etc.)

To analyze other geometries you need to modify one of the input files.

