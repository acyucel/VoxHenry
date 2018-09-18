clc
close all
clear all

% use recent DIRECTFN (faster)
use_recent_DIRECTFN=1;

filesep = filesep();

% compile MEX-based matrix-fill routines
disp('Compiling matrix-fill routines...')

if (use_recent_DIRECTFN)
  cd(['DIRECTFN', filesep, 'mex']);
  build
else
  if (ispc)
    cd src_lin_vie\singular\singular_win_compiled\
  else
    cd src_lin_vie/singular/singular_linux_compiled/
  end

  % if running under Octave, use the mex build specific linear_build
  % remark: not ported for Linux yet!
  if(exist ('OCTAVE_VERSION', 'builtin') > 0)
      linear_build_octave
  else
      linear_build
  end
end

cd ..
cd ..
disp('Done... Compiling matrix-fill routines')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('******************************************************************')
disp('Attention: VoxHenry uses ldlt decomposition of CHOLMOD in suitesparse,')
disp('which requires minimum memory for the Schur complement of sparse preconditioner')
disp('while providing minimum CPU times for multiplication and decomposition!')
disp('Although compiled CHOLMOD modules of suitesparse is provided in ')
disp('suitesparse_cholmod_XXX_compiled folders, where XXX is win or linux,')
disp('users are encouraged to compile the latest suitesparse, which can be')
disp('downloaded from http://faculty.cse.tamu.edu/davis/suitesparse.html')
disp('If suitesparse CANNOT be used in the current configuration, set "fl_cholmod"')
disp('flag in lse_sparse_precon_prepare to 0, for using ldlt decomposition of Matlab,')
disp('even though it is not as optimal as that of CHOLMOD. ')
disp('Alternatively, users can use lu decomposition of Matlab instead ldlt, which ')
disp('requires more memory, but runs as fast as ldlt of CHOLMOD. To do that, users should set the')
disp('slct_decomp_sch flag in lse_sparse_precon_prepare to lu_decomp.')
disp('******************************************************************')
% install Suitesparse for sparse preconditioner

% disp('Compiling Suitesparse routines...')
% cd suitesparse_cholmod\
% cholmod_install
% cd ..