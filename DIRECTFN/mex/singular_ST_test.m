function I_ST = singular_ST_test(volume_ker,l,k0,dx,Np,ker_type,np,nq,r1,r2,r3,r4)

% function is called with l-1 due to the difference in indexing between
% matlab and C++.

 k0 = 0.2*pi;
 Np = 5;
 dx = 1.0;

 % test 1
 r1 = [0., 0., 0.];
 r2 = [0.,dx,0.];
 r3 = [0., dx, dx];
 r4 = [0., 0., dx];
 nq = [1.,0.,0.];
 np = [1.,0.,0.];


% centre of squares
rq_c = (r1 + r3)/2;
rp_c = (r1 + r3)/2;


for volume_ker = 2:3
  for ker_type = 1:4
    for l = 1:9
I_ST = directfn_quad_st_plan_voxhenry(r1,r2,r3,r4,Np,Np,Np,Np,k0,dx,rq_c,rq_c,nq,np,ker_type,l,volume_ker)
    end
  end
end


end