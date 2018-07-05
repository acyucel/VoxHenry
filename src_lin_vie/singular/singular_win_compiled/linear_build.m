clc; close all; clear all;
if(ispc)
    str = 'mex ''-I..\include'' -c -outdir ''..\'' ';
else
    str = 'mex ''-I../include'' -c -outdir ''../'' ';
end

if(ispc)
    cd '.\utils';
else
    cd './utils';
end
list = strjoin(cellstr(ls('*.cpp')));
if(ispc)
    eval([str list]);
else
    %unix([str list]);
    mex('-c','-I../include','-outdir','../','*.cpp')
end

if(ispc)
    cd '..\WS_EA';
else
    cd '../WS_EA';
end
list = strjoin(cellstr(ls('*.cpp')));
if(ispc)
    eval([str list]);
else
    %unix([str list]);
    mex('-c','-I../include','-outdir','../','*.cpp')
end

if(ispc)
    cd '..\WS_ST';
else
    cd '../WS_ST';
end
list = strjoin(cellstr(ls('*.cpp')));
if(ispc)
    eval([str list]);
else
    %unix([str list]);
    mex('-c','-I../include','-outdir','../','*.cpp')
end

if(ispc)
    cd '..\WS_VA';
else
    cd '../WS_VA';
end
list = strjoin(cellstr(ls('*.cpp')));

if(ispc)
    eval([str list]);
else
    %unix([str list]);
    mex('-c','-I../include','-outdir','../','*.cpp')
end


if(ispc)
    cd '..';
else
    cd ..
end
 
if(ispc) 
    list = strjoin(cellstr(ls('*.obj')));
else
    list = strjoin(cellstr(ls('*.o')));
end

if(ispc)
    eval(['mex -output solve_ea ' list ' solve_ea_mex.cpp create_EA.cpp Kernels.cpp'])
    eval(['mex -output solve_va ' list ' solve_va_mex.cpp create_VA.cpp Kernels.cpp'])
    eval(['mex -output solve_st ' list ' solve_st_mex.cpp create_ST.cpp Kernels.cpp'])
else
    mex -output solve_ea *.o solve_ea_mex.cpp create_EA.cpp Kernels.cpp
    mex -output solve_va *.o solve_va_mex.cpp create_VA.cpp Kernels.cpp
    mex -output solve_st *.o solve_st_mex.cpp create_ST.cpp Kernels.cpp
    %unix(['mex -output solve_ea ' list ' solve_ea_mex.cpp create_EA.cpp Kernels.cpp'])
    %unix(['mex -output solve_va ' list ' solve_va_mex.cpp create_VA.cpp Kernels.cpp'])
    %unix(['mex -output solve_st ' list ' solve_st_mex.cpp create_ST.cpp Kernels.cpp'])
end

if(ispc) 
    delete('*.obj');
else
    delete('*.o');
end