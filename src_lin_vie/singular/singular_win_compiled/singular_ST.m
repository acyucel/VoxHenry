function I_ST = singular_ST(volume_ker,l,k0,dx,Np,ker_type,np,nq,r1,r2,r3,r4)

% function is called with l-1 due to the difference in indexing between
% matlab and C++.

% k0 = 0.2*pi;
% Np = 5;
% dx = 1;

% r1 = [0., 0., 0.];
% r2 = [0.,dx,0.];
% r3 = [0., dx, dx];
% r4 = [0., 0., dx];
% r5 = [0., 0., 2*dx];
% r6 = [0., dx, 2*dx];
% nq = [0.,0.,1.];
% np = [0.,0.,1.];
% rq_c = [dx/2,dx/2,0.];
% rp_c = [dx/2, 3*dx/2, 0.];

% centre of squares
rq_c = (r1 + r3)/2;

I_ST = solve_st(r1,r2,r3,r4,Np,Np,Np,Np,k0,dx,rq_c,rq_c,nq,np,ker_type,l,volume_ker);



end