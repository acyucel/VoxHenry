function I_ST = singular_ST(volume_ker,l,k0,dx,Np,ker_type,np,nq,r1,r2,r3,r4)

% centre of squares
% actually not needed, input parameters kept for compatibility
rq_c = (r1 + r3)/2;

I_ST = directfn_quad_st_plan_voxhenry(r1,r2,r3,r4,Np,Np,Np,Np,k0,dx,rq_c,rq_c,nq,np,ker_type,l,volume_ker);

% debug
%disp(['r1 ', num2str(r1'), ' r2 ', num2str(r2'), ' r3 ', num2str(r3'), ' r4 ', num2str(r4'), ' rq_c ', num2str(rq_c'), ' nq ', num2str(nq'), ' np ', num2str(np')])
%disp(['I_ST  using DIRECTFN ', num2str(I_ST)])

end
