addpath('../../mex')

d = 0.1;

k0 = 2*pi;


N1 = 10;
N2 = 10;
N3 = 10;
N4 = 10;

r1 = [ 0.0 , 0.0 , 0.0 ];
r2 = [ d , 0.0 , 0.0 ];
r3 = [ d,  d , 0.0 ];
r4 = [ 0.0 , d , 0.0 ];
r5 = [ 0.0, 2 * d, 0.0 ];
r6 = [ d, 2 * d, 0.0 ];

I_plan_const = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, N1, N2, N3, N4, k0, 'Constant') 

I_plan_WS = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, N1, N2, N3, N4, k0, 'Vector_WS') 

I_plan_SS = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, N1, N2, N3, N4, k0, 'Vector_SS') 



dd = d / sqrt(2.0);

r1 = [ d, 0.0, 0.0 ];
r2 = [ d, 0.0, 0.5*d ];
r3 = [ d, 0.0, d ];
r4 = [ dd, dd, 0.0 ];
r5 = [ dd, dd, 0.5*d ];
r6 = [ dd, dd, d ];
r7 = [ 0.0, d, 0.0 ];
r8 = [ 0.0, d, 0.5*d ];
r9 = [ 0.0, d, d ];
r10 = [ -dd, dd, 0.0 ];
r11 = [ -dd, dd, 0.5*d ];
r12 = [ -dd, dd, d ];
r13 = [ -d, 0.0, 0.0 ];
r14 = [ -d, 0.0, 0.5*d ];
r15 = [ -d, 0.0, d ];

I_curv_const = directfn_quad_ea_curv(r1, r2, r3, r4, r5, r6, ...
    r7, r8, r9 , r10, r11, r12, r13, r14, r15, N1, N2, N3, N4, k0, 'Constant') 

I_curv_WS = directfn_quad_ea_curv(r1, r2, r3, r4, r5, r6, ...
    r7, r8, r9, r10, r11, r12, r13, r14, r15,N1, N2, N3, N4, k0, 'Vector_WS') 


I_curv_SS = directfn_quad_ea_curv(r1, r2, r3, r4, r5, r6, ...
    r7, r8, r9, r10, r11, r12, r13, r14, r15, N1, N2, N3, N4, k0, 'Vector_SS') 
