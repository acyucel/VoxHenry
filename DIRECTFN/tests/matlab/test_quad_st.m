addpath('../../mex')

d = 0.1;

k0 = 2*pi;


N1 = 10;
N2 = 10;
N3 = 10;
N4 = 10;


r1 = [0.0 , 0.4*d , 0.0 ];
r2 = [ d , 0.0 , 0.0 ];
r3 = [0.8*d,  d , 0.0 ];
r4 = [ 0.2*d , 0.8*d , 0.0 ];

I_plan_const = directfn_quad_st_plan(r1, r2, r3, r4, N1, N2, N3, N4, k0, 'Constant') 

I_plan_WS = directfn_quad_st_plan(r1, r2, r3, r4, N1, N2, N3, N4, k0, 'Vector_WS') 

dd = d / sqrt(2.0);

r1 = [ 0.0, d, 0.0 ];
r2 = [ d,0.0,0.0 ];
r3 =[ d,0.0,d ];
r4 =[ 0.0, d, d ];
r5 =[ dd,dd,0.0 ];
r6 =[ d, 0.0, 0.5*d ];
r7 =[ dd, dd, d ];
r8 = [ 0.0,d,0.5*d ];
r9 =[ dd,dd,0.5*d ];

I_curv_const = directfn_quad_st_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, N1, N2, N3, N4, k0, 'Constant')
I_curv_WS = directfn_quad_st_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, N1, N2, N3, N4, k0, 'Vector_WS')