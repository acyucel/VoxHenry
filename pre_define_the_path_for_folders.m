% use recent DIRECTFN (faster)
use_recent_DIRECTFN=1;

% find the current folder
currentfolder = pwd;

% get system-dependent file separator
filesep = filesep();
% get system-dependent path separator
psep = pathsep();

% include all the needed path (addpath() adds recursively),
% but do NOT include the Git directory

if(use_recent_DIRECTFN)
  p1  = genpath([currentfolder, filesep, 'DIRECTFN']);
  addpath(p1);
end

p2  = genpath([currentfolder, filesep, 'Results']);
p3  = genpath([currentfolder, filesep, 'results_numex1_straight_conductor']);
p4  = genpath([currentfolder, filesep, 'results_numex2_wire']);
p5  = genpath([currentfolder, filesep, 'results_numex3_circular_coil']);
p6  = genpath([currentfolder, filesep, 'results_numex4_square_coil']);
p7  = genpath([currentfolder, filesep, 'src_iterative_solver']);
p8  = genpath([currentfolder, filesep, 'src_lin_vie']);
p9  = genpath([currentfolder, filesep, 'src_lse_formation']);
p10 = genpath([currentfolder, filesep, 'src_post_process']);
p11 = genpath([currentfolder, filesep, 'src_pre_process']);
addpath([p2, psep, p3, psep, p4, psep, p5, psep, p6, psep, p7, psep, p8, psep, p9, psep, p10, psep, p11]);

if (ispc)
  p12 = genpath([currentfolder, filesep, 'suitesparse_cholmod_win_compiled']);
else
  p12 = genpath([currentfolder, filesep, 'suitesparse_cholmod_linux_compiled']);
end
addpath(p12);

if(use_recent_DIRECTFN)
  rmpath([currentfolder, filesep, 'src_lin_vie', filesep, 'singular', filesep, 'singular_win_compiled']);
  rmpath([currentfolder, filesep, 'src_lin_vie', filesep, 'singular', filesep, 'singular_linux_compiled']);
else
  if (ispc) 
    rmpath([currentfolder, filesep, 'src_lin_vie', filesep, 'singular', filesep, 'singular_linux_compiled']);
  else
    rmpath([currentfolder, filesep, 'src_lin_vie', filesep, 'singular', filesep, 'singular_win_compiled']);
  end
end

  

    
  
