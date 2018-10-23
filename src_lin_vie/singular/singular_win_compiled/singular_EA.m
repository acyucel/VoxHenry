function I_EA = singular_EA(volume_ker,l,k0,dx,Np,ker_type,np,nq,r1,r2,r3,r4,r5,r6)

% centre of squares
rq_c = (r1 + r3)/2;
rp_c = (r4 + r6)/2;

I_EA = solve_ea(r1,r2,r3,r4,r5,r6,Np,Np,Np,Np,k0,dx,rq_c,rp_c,nq,np,ker_type,l,volume_ker);

% debug
%disp(['I_EA ', num2str(I_EA)])

end