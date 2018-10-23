
lib_path = fullfile('-L..','lib','unix');
%lib_path = fullfile('-L..','lib','win');
%lib_path = fullfile('-L..','lib');
inc_path = fullfile('-I..','include');
addpath('../')

eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_quad_st_plan.cpp']);
eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_quad_st_curv.cpp']);
eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_quad_ea_plan.cpp']);
eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_quad_ea_curv.cpp']);
eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_quad_va_plan.cpp']);
eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_quad_va_curv.cpp']);

eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_tri_st_plan.cpp']);
eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_tri_ea_plan.cpp']);
eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_tri_va_plan.cpp']);

eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_quad_st_plan_voxhenry.cpp']);
eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_quad_ea_plan_voxhenry.cpp']);
eval(['mex -v ' ' ' inc_path ' ' lib_path ' ' '-ldirectfn' ' ' 'directfn_quad_va_plan_voxhenry.cpp']);

