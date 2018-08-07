% use recent DIRECTFN (faster)
use_recent_DIRECTFN=1;

% find the current folder
currentfolder = pwd;
teracisfolder=currentfolder;

% obtain the string with the recursive paths, but do NOT include the Git directory
p = genpath(teracisfolder, '.git');
   
% add the path tree to the current path
addpath(p);

% get system-dependent file separator
filesep = filesep();

currentfolder = pwd;

if(use_recent_DIRECTFN)
  rmpath([currentfolder, filesep, 'src_lin_vie', filesep, 'singular', filesep, 'singular_win_compiled']);
  rmpath([currentfolder, filesep, 'src_lin_vie', filesep, 'singular', filesep, 'singular_linux_compiled']);
else
  rmpath([currentfolder, filesep, 'DIRECTFN', filesep, 'mex']);
  if (ispc) 
    rmpath([currentfolder, filesep, 'src_lin_vie', filesep, 'singular', filesep, 'singular_linux_compiled']);
  else
    rmpath([currentfolder, filesep, 'src_lin_vie', filesep, 'singular', filesep, 'singular_win_compiled']);
  end
end

if (ispc)
  % remove linux suitsparse
  rmpath([currentfolder, filesep, 'suitesparse_cholmod_linux_compiled']);
else
  rmpath([currentfolder, filesep, 'suitesparse_cholmod_win_compiled']);
end
  

    
  
