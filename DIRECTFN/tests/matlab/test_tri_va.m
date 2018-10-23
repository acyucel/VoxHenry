addpath('../../mex')

d = 0.1;

k0 = 2*pi;
r1 = [ 0.0, d, 0.0 ];
r2 = [ 0.0, 0.0, 0.0 ];
r3 = [ d, 0.0, 0.0 ];
r4 = [ d, d, d ];
r5 = [ 2*d, d, d ];


N1 = 10;
N2 = 10;
N3 = 10;
N4 = 10;

I_const = directfn_tri_va_plan(r1, r2, r3, r4, r5, N1, N2, N3, N4, k0, 'Constant')
I_RWG_WS = directfn_tri_va_plan(r1, r2, r3, r4, r5, N1, N2, N3, N4, k0, 'RWG_WS')
I_RWG_SS = directfn_tri_va_plan(r1, r2, r3, r4, r5, N1, N2, N3, N4, k0, 'RWG_SS')
I_nxRWG_SS = directfn_tri_va_plan(r1, r2, r3, r4, r5, N1, N2, N3, N4, k0, 'nxRWG_SS')



