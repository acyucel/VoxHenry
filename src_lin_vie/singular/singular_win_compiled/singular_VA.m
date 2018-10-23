function I_VA = singular_VA(volume_ker,l,k0,dx,Np,ker_type,np,nq,r1,r2,r3,r4,r5,r6,r7)

% centre of squares
rq_c = (r1 + r3)/2;
rp_c = (r3 + r6)/2;

I_VA = solve_va(r1,r2,r3,r4,r5,r6,r7,Np,Np,Np,Np,k0,dx,rq_c,rp_c,nq,np,ker_type,l,volume_ker);

% debug
%disp(['I_VA ', num2str(I_VA)])

end