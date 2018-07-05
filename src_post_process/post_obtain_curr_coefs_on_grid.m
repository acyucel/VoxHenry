function [Jx_currs_grid,Jy_currs_grid,Jz_currs_grid,J2d_currs_grid,J3d_currs_grid,cmin,cmax]=post_obtain_curr_coefs_on_grid(x,Mc)

[L,M,N] = size(Mc); % domain size

tola=1e-12;
dum=1;
ignore_unk=zeros(L*M*N,1);
for mz = 1:N;
    for my = 1:M;
        for mx = 1:L;
            if(abs(Mc(mx,my,mz)) < tola)
                ignore_unk(dum) = 1;
            end
            dum=dum+1;
        end;
    end;
end

num_cube=L*M*N;
num_cube_new=num_cube-sum(ignore_unk);

Jx_currs=x(1:num_cube_new);Jy_currs=x(num_cube_new+1:2*num_cube_new); Jz_currs=x(2*num_cube_new+1:3*num_cube_new);
J2d_currs=x(3*num_cube_new+1:4*num_cube_new);J3d_currs=x(4*num_cube_new+1:5*num_cube_new);

cmin=min(abs([Jx_currs;Jy_currs;Jz_currs;]));
cmax=max(abs([Jx_currs;Jy_currs;Jz_currs;]));

Jx_currs_grid=zeros(L,M,N);Jy_currs_grid=zeros(L,M,N); Jz_currs_grid=zeros(L,M,N);
J2d_currs_grid=zeros(L,M,N);J3d_currs_grid=zeros(L,M,N);

dum=1;dum2=1;
for mz = 1:N
    for my = 1:M
        for mx = 1:L
            if ( ignore_unk(dum) == 0 );
               Jx_currs_grid(mx,my,mz)= Jx_currs(dum2);
               Jy_currs_grid(mx,my,mz)= Jy_currs(dum2);
               Jz_currs_grid(mx,my,mz)= Jz_currs(dum2);
               J2d_currs_grid(mx,my,mz)= J2d_currs(dum2);
               J3d_currs_grid(mx,my,mz)= J3d_currs(dum2);
               dum2=dum2+1;
            end
            dum=dum+1;
        end
    end
end