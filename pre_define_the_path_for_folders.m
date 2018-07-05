% find the current folder
currentfolder = pwd;
teracisfolder=currentfolder;

% obtain the string with the recursive paths
p = genpath(teracisfolder);

% add the path tree to the current path
addpath(p);

if (ispc)
    currentfolder = pwd;
    rmpath([currentfolder,'\suitesparse_cholmod_linux_compiled']);
    rmpath([currentfolder,'\src_lin_vie\singular\singular_linux_compiled']);
else
    currentfolder = pwd;
    rmpath([currentfolder,'/suitesparse_cholmod_win_compiled']);
    rmpath([currentfolder,'/src_lin_vie/singular/singular_win_compiled']);
end