function I_ST = singular_ST(volume_ker,l,k0,dx,Np,ker_type,np,nq,r1,r2,r3,r4)

% centre of squares
rq_c = (r1 + r3)/2;

I_ST = solve_st(r1,r2,r3,r4,Np,Np,Np,Np,k0,dx,rq_c,rq_c,nq,np,ker_type,l,volume_ker);

% debug
%disp(['I_ST ', num2str(I_ST)])

end