function I_EA = singular_EA(volume_ker,l,k0,dx,Np,ker_type,np,nq,r1,r2,r3,r4,r5,r6)

% volume_ker: the initial kernel of the volume volume integral G, 1/R, R
% that takes values 1,2,3 respectively 
% l: indexing for the integrals according to eq. 1 
% ker_type: takes values 1-4 referring to the each of the 4 surface-surface
% reduced kernels

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
rp_c = (r4 + r6)/2;

% call the mexed singular integral routines for edge adjacent case
I_EA = solve_ea(r1,r2,r3,r4,r5,r6,Np,Np,Np,Np,k0,dx,rq_c,rp_c,nq,np,ker_type,l,volume_ker);

end