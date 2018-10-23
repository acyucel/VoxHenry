function [fN_all,st_sparse_precon] = lse_generate_circulant_tensor(dx,ko,L,M,N,fl_no_fft)

% Code addition for running subroutine seperately
% clc; close all; clear all
%ko = 1;
%dx = 1e-6;
%L=3; M=3; N=1;
      
       
% if 'lse_enhanced_near' == 3, use pre-calculated near interactions when 
%     calculating near G_mn terms
% if 'lse_enhanced_near' == 2, exploit symmetry when calculating near G_mn terms
% if 'lse_enhanced_near' == 1, unroll the for loops to improve
%     parallelization for near interactions calculation
% if 'lse_enhanced_near' == 0, use the original brute force way to compute G_mn terms
lse_enhanced_near = 3;
% if 'lse_exploit_symmetry' == 0, use the original brute force way to compute G_mn terms
% if 'lse_exploit_symmetry' == 1, exploit symmetry when calculating G_mn terms
lse_exploit_symmetry = 1;

fl_oneoverr_kernel=1;
if(fl_oneoverr_kernel == 1)
    ko_grfn=6.287535065855046e-08; % ko of 3Hz
else
    ko_grfn=ko;
end


% 2) Computing interactions with Ioannis code
% 2a) Specify the number of integration points - far, near

% number of 1D points for far volume-volume integrals
% (the total points in the 6D quadrature will be Np_1D_far_V^6)
Np_1D_far_V    = 3;
% number of 1D points for medium volume-volume integrals
Np_1D_medium_V = 5;
% number of 1D points for near surface-surface and singular integrals 
% you can try to reduce the number of points to around 10, however for some
% kernels 12 points are required for convergence
Np_1D_near_S   = 12;

% Set region of medium distances cells for higher order quadrature
n_medium = 5; 

% centre of the source voxel - reference cell
r_n = [0 0 0]';
% distance vector
dy = dx; dz = dx; 
d = [dx,dy,dz];

% 2b) Compute interaction between far distance cells - the volume-volume integrals
G_mn = zeros(L,M,N,7);
infofN = whos('G_mn'); memestimated = 2*infofN.bytes/(1024*1024);
disp(['  Memory for temporarily storing G_mn (MB) ::: ' , num2str(memestimated)]);

Np = Np_1D_far_V;
% 6D Clenshaw - Curtis weights and points
[W6D,X,Y,Z,Xp,Yp,Zp] = weights_points(Np,6);
disp(['  Start computing far interactions'])
tic
parfor mx = 1:L
    for my = 1:M
        for mz = 1:N

            m = [mx,my,mz];
            % centre of the observation voxel
            r_m = ((m-1) .* d)';

            %G_mn(mx,my,mz,:) = volume_volume(W6D,X,Y,Z,Xp,Yp,Zp,ko_grfn,r_m,r_n,dx);
            G_mn(mx,my,mz,:) = volume_volume(W6D,X,Y,Z,Xp,Yp,Zp,ko_grfn,r_m,r_n,dx,1);

        end
    end
end
disp(['  Time for computing far interactions ::: ',num2str(toc)])

% 2d) Compute interaction between medium distance cells - the surface-surface integrals

Np = Np_1D_medium_V;
% 6D Clenshaw - Curtis weights and points
[W6D,X,Y,Z,Xp,Yp,Zp] = weights_points(Np,6);

n_mediumL = n_medium;
n_mediumM = n_medium;
n_mediumN = n_medium;
% adjust domain for singular integrals calculations
if L <n_mediumL; n_mediumL = L; end
if M <n_mediumM; n_mediumM = M; end
if N <n_mediumN; n_mediumN = N; end

% if should not exploit symmetry, or if not a cubic volume, use standard method
if(lse_exploit_symmetry == 0 || n_mediumL ~= n_mediumM || n_mediumL ~= n_mediumN ) 
    tic
    parfor mx = 1:n_mediumL
        for my = 1:n_mediumM
            for mz = 1:n_mediumN
                
                m = [mx,my,mz];
                % centre of the observation voxel
                r_m = ((m-1) .* d)';
                
                G_mn(mx,my,mz,:) = volume_volume(W6D,X,Y,Z,Xp,Yp,Zp,ko_grfn,r_m,r_n,dx,1);
                
            end
        end
    end
    disp(['  Time for computing medium interactions ::: ',num2str(toc)])
else
    % compute medium interactions exploiting symmetry
    % remark: volume must be cubic, i.e. 'mediumL' = 'mediumM' = 'mediumN'
    % remark: MUST zero the elements 6,9,10 in the last mode of the tensor,
    %         as these calculations use accumulation 
    %         (e.g. G_mn(mx,my,mz,6) = G_mn(mx,my,mz,6) + K(3);) 
    %
    % Explanation of the symmetry exploit:
    % 
    % We observe that the voxel region is actually 1/8 of the 360 degrees
    % solid angle around the axis origin, and that calculating the 
    % integrals we have either pulse functions or linear functions
    % along x, y or z.
    % Then, rotating the coordinate reference, what was calculated for 
    % linear functions along x is now valid for y:
    %      z   x                                y   z
    %      |  /       becomes after rotation    |  /
    %      | /                                  | /
    %      |/____ y                             |/____ x
    %
    % The same is also true for z with another rotation.
    % We need therefore to calculate the integrals only for the linear
    % functions along x.
    % 
    % We observe then that in the 1/8 volume considered there is also
    % symmetry w.r.t the coordinate x, in particular we can swap z and y,
    % and the result is the same:
    %      z   x                                y   x
    %      |  /       has a symmetry along x    |  /
    %      | /                                  | /
    %      |/____ y                             |/____ z
    %
    % So we only need to calculate integrals for x on about half of the 
    % voxels (the diagonal is included, so slightly more than the half)
    %
    % In the end, we also note from the integrals calculations 
    % that Gx,2D = -G2D,x so this saves yet another integral
     
    tic
    
    % let's get the indexes of the voxels that we need to consider
    % for the interactions. Note that we consider only half of the voxels 
    % (plus the diagonal)
    inter_idx = zeros( (n_mediumL * (n_mediumM + 1) * n_mediumN) / 2,1);
    dum = 1;
    for mx = 1:n_mediumL
        for my = 1:n_mediumM
            % to further exploit symmetry, let's calculate only helf of the quadrant
            for mz = 1:my
                inter_idx(dum) = sub2ind(size(G_mn), mx, my, mz);
                dum = dum + 1;
            end
        end
    end

    % Now let's create a temporary storage for the kernels calculated for
    % each voxel (three).
    % This is needed as MatLab does not like the (apparently) random 
    % insertion of elements in a matrix, plus some elements of G_mn
    % calculated via symmetry needs accumulation, so different threads
    % may try to access the same element at the same time, must avoid that
    K_tmp = zeros(size(inter_idx,1), 3);
    parfor dum = 1:size(inter_idx,1)
      
        [mx,my,mz] = ind2sub(size(G_mn),inter_idx(dum));
        m = [mx,my,mz];
        % centre of the observation voxel
        r_m = ((m-1) .* d)';

        K_tmp(dum,:) = volume_volume_sym(W6D,X,Y,Z,Xp,Yp,Zp,ko_grfn,r_m,r_n,dx,1,1);
    end
                
    % zero all elements within the medium interaction range (see remark above)
    G_mn(1:n_mediumL,1:n_mediumM,1:n_mediumN,:) = zeros(n_mediumL,n_mediumM,n_mediumN,7);
    % and now properly insert the near G_mn elements in the 'G_mn' tensor. 
    for dum = 1:size(inter_idx,1)
        [mx,my,mz] = ind2sub(size(G_mn),inter_idx(dum));
 
        % this is only half (1/8 of complete angle around the x-axis);
        % must copy to the other half, to cover the full quadrant
        
        % everything along x, first half of the quadrant
        
        % Gx,x
        G_mn(mx,my,mz,1) = K_tmp(dum,1);
        % G2D,x
        G_mn(mx,my,mz,2) = -K_tmp(dum,2);
        % G2D,2D
        % linear(x)-linear(x') part
        G_mn(mx,my,mz,4) = G_mn(mx,my,mz,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(x)-linear(x') part
        G_mn(mx,my,mz,6) = G_mn(mx,my,mz,6) + K_tmp(dum,3);
        % G3D,3D
        % linear(x)-linear(x') part
        G_mn(mx,my,mz,7) = G_mn(mx,my,mz,7) + K_tmp(dum,3);

        % everything along x, second half of the quadrant (exchanging y and z)
        
        % Gx,x
        G_mn(mx,mz,my,1) = K_tmp(dum,1);
        % G2D,x
        G_mn(mx,mz,my,2) = -K_tmp(dum,2);
        % G2D,2D
        % linear(x)-linear(x') part
        G_mn(mx,mz,my,4) = G_mn(mx,mz,my,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(x)-linear(x') part
        G_mn(mx,mz,my,6) = G_mn(mx,mz,my,6) + K_tmp(dum,3);
        % G3D,3D
        % linear(x)-linear(x') part
        G_mn(mx,mz,my,7) = G_mn(mx,mz,my,7) + K_tmp(dum,3);

        
        % everything along y

        % G2D,y
        G_mn(mz,mx,my,3) = K_tmp(dum,2);
        % G2D,2D
        % linear(y)-linear(y') part
        G_mn(mz,mx,my,4) = G_mn(mz,mx,my,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(y)-linear(y') part
        G_mn(mz,mx,my,6) = G_mn(mz,mx,my,6) - K_tmp(dum,3);
        % G3D,3D
        % linear(y)-linear(y') part
        G_mn(mz,mx,my,7) = G_mn(mz,mx,my,7) + K_tmp(dum,3);

        % everything along y, second half of the quadrant (exchanging x and z)

        % G2D,y
        G_mn(my,mx,mz,3) = K_tmp(dum,2);
        % G2D,2D
        % linear(y)-linear(y') part
        G_mn(my,mx,mz,4) = G_mn(my,mx,mz,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(y)-linear(y') part
        G_mn(my,mx,mz,6) = G_mn(my,mx,mz,6) - K_tmp(dum,3);
        % G3D,3D
        % linear(y)-linear(y') part
        G_mn(my,mx,mz,7) = G_mn(my,mx,mz,7) + K_tmp(dum,3);

        
        % everything along z

        % G3D,z
        G_mn(my,mz,mx,5) = 2*K_tmp(dum,2);
        % G3D,3D
        % linear(z)-linear(z') part
        G_mn(my,mz,mx,7) = G_mn(my,mz,mx,7) + 4*K_tmp(dum,3);

        % everything along z, second half of the quadrant (exchanging x and y)

        % G3D,z
        G_mn(mz,my,mx,5) = 2*K_tmp(dum,2);
        % G3D,3D
        % linear(z)-linear(z') part
        G_mn(mz,my,mx,7) = G_mn(mz,my,mx,7) + 4*K_tmp(dum,3);
        
    end
    disp(['  Time for computing medium interactions exploiting symmetry ::: ',num2str(toc)])
    % line commented out to use G_mn directly; even if in Matlab/Octave you might save
    % some time to use an intermediate 'G_mn_test' smaller tensor and copy it to 'G_mn' at the end
    %G_mn(1:n_mediumL,1:n_mediumL,1:n_mediumL,:) = G_mn_test(1:n_mediumL,1:n_mediumL,1:n_mediumL,:);
end

% 2e) Compute interaction between near distance cells - the surface-surface integrals

Np = Np_1D_near_S;
% 4D Gauss Legendre points and weights
[W4D,A,B,C,D] = weights_points(Np,4);

Ls = 2;Ms = 2;Ns = 2;
% adjust domain for singular integrals calculations
if L <2; Ls = 1; end
if M <2; Ms = 1; end
if N <2; Ns = 1;end

% calculate the constant coefficinets in front of the surface-surface
% integrals
[I1_co, I2_co, I3_co, I4_co]  = surface_surface_coeff(dx,ko_grfn);

if (lse_enhanced_near == 3 && Ls == Ms && Ls == Ns)
  
    % compute interactions using pre-calculated near interactions
    % and exploiting symmetry
    % remark: volume must be cubic, i.e. 'Ls' = 'Ms' = 'Ns'
    % remark: MUST zero the elements 6,9,10 in the last mode of the tensor,
    %         as these calculations use accumulation 
    %         (e.g. G_mn(mx,my,mz,6) = G_mn(mx,my,mz,6) + K(3);) 

    tic
    
    % let's get the indexes of the voxels that we need to consider
    % for the interactions. Note that we consider only half of the voxels 
    % (plus the diagonal). In 2D this is n*(n+1)/2, in 3D n*n*(n+1)/2 of course
    inter_idx = zeros( (Ls * (Ms + 1) * Ns) / 2,1);
    dum = 1;
    for mx = 1:Ls
        for my = 1:Ms
            % to further exploit symmetry, let's calculate only helf of the quadrant
            for mz = 1:my
                inter_idx(dum) = sub2ind(size(G_mn), mx, my, mz);
                dum = dum + 1;
            end
        end
    end

    % Now let's create a temporary storage for the kernels calculated for
    % each voxel (three).
    % This is needed as MatLab does not like the (apparently) random 
    % insertion of elements in a matrix, plus some elements of G_mn
    % calculated via symmetry needs accumulation, so different threads
    % may try to access the same element at the same time, must avoid that
    K_tmp = zeros(size(inter_idx,1), 3);

    % load the pre-calculated near interactions
    K_tmp = dlmread(['src_lin_vie', filesep(), 'K_near_interactions_1mm.txt']);
    
    if(size(inter_idx,1) ~= size(K_tmp,1))
        disp(['ERROR: loaded K-interactions matrix with length  ', num2str(size(K_load,1)), ' while we expected ', num2str(size(inter_idx,1))]);
    end
    
    % scale the K_tmp terms. The method leverages the fact that the kernel has
    % the special property that g(k*r) = f(k) * g(r).
    % In 1D, we use this fact to transform the integral 
    %
    %  /kb
    %  | g(x)dx 
    % /ka
    %
    % to a standard form
    %
    %  /a
    %  | g(y)dy
    % /b
    %
    % So transforming x = k*y we have
    %
    %  /kb        /b               /b                         /b
    %  | g(x)dx = | g(k*y)*k*dy =  | f(k)*g(y)*k*dy = f(k)*k* | g(y)dy
    % /ka        /a               /a                         /a  
    %
    % In our particular 3D case, and with the kernel 1/r, we have g(k*r) = 1/k * g(r)
    % while the Jacobian is k^3; plus, we have a double integral (two Jacobians),
    % that gives a k^3 * k^3 factor, therefore we multiply by k^3*k^3/k = k^5
    % The same works for the kernel x*1/r, where g(k*r) = k/k * g(r) = g(r)
    % so we multiply by k^6
    % And of course for the kernel x*x'*1/r, where g(r*k) = k^2/k * g(r) = k*g(r)
    % so we multiply by k^7
    
    K_scale = dx / 1e-3;
    K_tmp(:,1) = K_tmp(:,1) * K_scale ^ 5;
    K_tmp(:,2) = K_tmp(:,2) * K_scale ^ 6;
    K_tmp(:,3) = K_tmp(:,3) * K_scale ^ 7;
     
    % zero all elements within the near interaction range (see remark above)
    G_mn(1:Ls,1:Ms,1:Ns,:) = zeros(Ls,Ms,Ns,7);
    % and now properly insert the near G_mn elements in the 'G_mn' tensor. 
    for dum = 1:size(inter_idx,1)
        [mx,my,mz] = ind2sub(size(G_mn),inter_idx(dum));
 
        % this is only half (1/8 of complete angle around the x-axis);
        % must copy to the other half, to cover the full quadrant
        
        % everything along x, first half of the quadrant
        
        % Gx,x
        G_mn(mx,my,mz,1) = K_tmp(dum,1);
        % G2D,x
        G_mn(mx,my,mz,2) = -K_tmp(dum,2);
        % G2D,2D
        % linear(x)-linear(x') part
        G_mn(mx,my,mz,4) = G_mn(mx,my,mz,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(x)-linear(x') part
        G_mn(mx,my,mz,6) = G_mn(mx,my,mz,6) + K_tmp(dum,3);
        % G3D,3D
        % linear(x)-linear(x') part
        G_mn(mx,my,mz,7) = G_mn(mx,my,mz,7) + K_tmp(dum,3);

        % everything along x, second half of the quadrant (exchanging y and z)
        
        % Gx,x
        G_mn(mx,mz,my,1) = K_tmp(dum,1);
        % G2D,x
        G_mn(mx,mz,my,2) = -K_tmp(dum,2);
        % G2D,2D
        % linear(x)-linear(x') part
        G_mn(mx,mz,my,4) = G_mn(mx,mz,my,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(x)-linear(x') part
        G_mn(mx,mz,my,6) = G_mn(mx,mz,my,6) + K_tmp(dum,3);
        % G3D,3D
        % linear(x)-linear(x') part
        G_mn(mx,mz,my,7) = G_mn(mx,mz,my,7) + K_tmp(dum,3);

        
        % everything along y

        % G2D,y
        G_mn(mz,mx,my,3) = K_tmp(dum,2);
        % G2D,2D
        % linear(y)-linear(y') part
        G_mn(mz,mx,my,4) = G_mn(mz,mx,my,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(y)-linear(y') part
        G_mn(mz,mx,my,6) = G_mn(mz,mx,my,6) - K_tmp(dum,3);
        % G3D,3D
        % linear(y)-linear(y') part
        G_mn(mz,mx,my,7) = G_mn(mz,mx,my,7) + K_tmp(dum,3);

        % everything along y, second half of the quadrant (exchanging x and z)

        % G2D,y
        G_mn(my,mx,mz,3) = K_tmp(dum,2);
        % G2D,2D
        % linear(y)-linear(y') part
        G_mn(my,mx,mz,4) = G_mn(my,mx,mz,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(y)-linear(y') part
        G_mn(my,mx,mz,6) = G_mn(my,mx,mz,6) - K_tmp(dum,3);
        % G3D,3D
        % linear(y)-linear(y') part
        G_mn(my,mx,mz,7) = G_mn(my,mx,mz,7) + K_tmp(dum,3);

        
        % everything along z

        % G3D,z
        G_mn(my,mz,mx,5) = 2*K_tmp(dum,2);
        % G3D,3D
        % linear(z)-linear(z') part
        G_mn(my,mz,mx,7) = G_mn(my,mz,mx,7) + 4*K_tmp(dum,3);

        % everything along z, second half of the quadrant (exchanging x and y)

        % G3D,z
        G_mn(mz,my,mx,5) = 2*K_tmp(dum,2);
        % G3D,3D
        % linear(z)-linear(z') part
        G_mn(mz,my,mx,7) = G_mn(mz,my,mx,7) + 4*K_tmp(dum,3);
        
    end
    
    disp(['  Time for computing near interactions ::: ',num2str(toc)])
    
elseif (lse_enhanced_near == 2 && Ls == Ms && Ls == Ns)
  
    % compute interactions exploiting symmetry
    % remark: volume must be cubic, i.e. 'Ls' = 'Ms' = 'Ns'
    % remark: MUST zero the elements 6,9,10 in the last mode of the tensor,
    %         as these calculations use accumulation 
    %         (e.g. G_mn(mx,my,mz,6) = G_mn(mx,my,mz,6) + K(3);) 

    tic
    
    % let's get the indexes of the voxels that we need to consider
    % for the interactions. Note that we consider only half of the voxels 
    % (plus the diagonal). In 2D this is n*(n+1)/2, in 3D n*n*(n+1)/2 of course
    inter_idx = zeros( (Ls * (Ms + 1) * Ns) / 2,1);
    dum = 1;
    for mx = 1:Ls
        for my = 1:Ms
            % to further exploit symmetry, let's calculate only helf of the quadrant
            for mz = 1:my
                inter_idx(dum) = sub2ind(size(G_mn), mx, my, mz);
                dum = dum + 1;
            end
        end
    end

    % Now let's create a temporary storage for the kernels calculated for
    % each voxel (three).
    % This is needed as MatLab does not like the (apparently) random 
    % insertion of elements in a matrix, plus some elements of G_mn
    % calculated via symmetry needs accumulation, so different threads
    % may try to access the same element at the same time, must avoid that
    K_tmp = zeros(size(inter_idx,1), 3);
    parfor dum = 1:size(inter_idx,1)
      
        [mx,my,mz] = ind2sub(size(G_mn),inter_idx(dum));
        m = [mx,my,mz];
        % centre of the observation voxel
        r_m = ((m-1) .* d)';

        % 1) calculate near interactions with singularity subtraction
        % method. VV for the smoothed kernel, SS for 1/R, SS for R/2
        I_V_smooth = volume_volume_sym(W6D,X,Y,Z,Xp,Yp,Zp,ko_grfn,r_m,r_n,dx,4,1);
        I_S2 = surface_surface_kernels_sym(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,2,1);
        I_S3 = surface_surface_kernels_sym(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,3,1);

        K_tmp(dum,:) = I_V_smooth + I_S2 + I_S3;
        
    end
    
    % export the 'K_tmp' array for re-using in fastest near interaction computation
    %dlmwrite(['src_lin_vie', filesep(), 'K_near_interactions_1mm.txt'], K_tmp, 'precision', 16);
                
    % zero all elements within the near interaction range (see remark above)
    G_mn(1:Ls,1:Ms,1:Ns,:) = zeros(Ls,Ms,Ns,7);
    % and now properly insert the near G_mn elements in the 'G_mn' tensor. 
    for dum = 1:size(inter_idx,1)
        [mx,my,mz] = ind2sub(size(G_mn),inter_idx(dum));
 
        % this is only half (1/8 of complete angle around the x-axis);
        % must copy to the other half, to cover the full quadrant
        
        % everything along x, first half of the quadrant
        
        % Gx,x
        G_mn(mx,my,mz,1) = K_tmp(dum,1);
        % G2D,x
        G_mn(mx,my,mz,2) = -K_tmp(dum,2);
        % G2D,2D
        % linear(x)-linear(x') part
        G_mn(mx,my,mz,4) = G_mn(mx,my,mz,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(x)-linear(x') part
        G_mn(mx,my,mz,6) = G_mn(mx,my,mz,6) + K_tmp(dum,3);
        % G3D,3D
        % linear(x)-linear(x') part
        G_mn(mx,my,mz,7) = G_mn(mx,my,mz,7) + K_tmp(dum,3);

        % everything along x, second half of the quadrant (exchanging y and z)
        
        % Gx,x
        G_mn(mx,mz,my,1) = K_tmp(dum,1);
        % G2D,x
        G_mn(mx,mz,my,2) = -K_tmp(dum,2);
        % G2D,2D
        % linear(x)-linear(x') part
        G_mn(mx,mz,my,4) = G_mn(mx,mz,my,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(x)-linear(x') part
        G_mn(mx,mz,my,6) = G_mn(mx,mz,my,6) + K_tmp(dum,3);
        % G3D,3D
        % linear(x)-linear(x') part
        G_mn(mx,mz,my,7) = G_mn(mx,mz,my,7) + K_tmp(dum,3);

        
        % everything along y

        % G2D,y
        G_mn(mz,mx,my,3) = K_tmp(dum,2);
        % G2D,2D
        % linear(y)-linear(y') part
        G_mn(mz,mx,my,4) = G_mn(mz,mx,my,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(y)-linear(y') part
        G_mn(mz,mx,my,6) = G_mn(mz,mx,my,6) - K_tmp(dum,3);
        % G3D,3D
        % linear(y)-linear(y') part
        G_mn(mz,mx,my,7) = G_mn(mz,mx,my,7) + K_tmp(dum,3);

        % everything along y, second half of the quadrant (exchanging x and z)

        % G2D,y
        G_mn(my,mx,mz,3) = K_tmp(dum,2);
        % G2D,2D
        % linear(y)-linear(y') part
        G_mn(my,mx,mz,4) = G_mn(my,mx,mz,4) + K_tmp(dum,3);
        % G2D,3D
        % linear(y)-linear(y') part
        G_mn(my,mx,mz,6) = G_mn(my,mx,mz,6) - K_tmp(dum,3);
        % G3D,3D
        % linear(y)-linear(y') part
        G_mn(my,mx,mz,7) = G_mn(my,mx,mz,7) + K_tmp(dum,3);

        
        % everything along z

        % G3D,z
        G_mn(my,mz,mx,5) = 2*K_tmp(dum,2);
        % G3D,3D
        % linear(z)-linear(z') part
        G_mn(my,mz,mx,7) = G_mn(my,mz,mx,7) + 4*K_tmp(dum,3);

        % everything along z, second half of the quadrant (exchanging x and y)

        % G3D,z
        G_mn(mz,my,mx,5) = 2*K_tmp(dum,2);
        % G3D,3D
        % linear(z)-linear(z') part
        G_mn(mz,my,mx,7) = G_mn(mz,my,mx,7) + 4*K_tmp(dum,3);
        
    end
    
    disp(['  Time for computing near interactions exploiting symmetry ::: ',num2str(toc)])
    
elseif (lse_enhanced_near == 1)
    % Parallelization requires the maximum number of threads working in
    % parallel. So if we use nested 'for' loops, only the external
    % one can be parallelized, and this is the ultimate limit of 
    % number of active threads. If the external for loop has less
    % iterations than the number of threads, there will be processors
    % sitting idle. The solution is to serialize the computation in
    % a single (parallel) for loop.
    % Note that the ultimate limit here is Ls * Ms * Ns number of
    % threads. If Ls, Ms, Ns are equal to 2, the limit is then 8 threads.
    % We can further enhance the parallelization for more threads
    % by splitting the job for calculating 'I_V_smooth', 'I_S2' and 'I_S3'
    % in parallel as well.

    % let's get the indexes of the voxels that we need to consider
    % for near interactions
    near_idx = zeros(Ls * Ms * Ns,1);
    dum = 1;
    for mx = 1:Ls
        for my = 1:Ms
            for mz = 1:Ns
                near_idx(dum) = sub2ind(size(G_mn), mx, my, mz);
                dum = dum + 1;
            end
        end
    end
    % now let's create a temporary storage for the G_mn near elements.
    % this is needed as MatLab does not like the (apparently) random 
    % insertion of elements in G_mn. We know that this is instead 
    % not an issue (no dependances between the threads),
    % but need to work it around here.
    G_mn_tmp = zeros(Ls * Ms * Ns, 7);
    parfor dum = 1:size(near_idx,1)

                [mx,my,mz] = ind2sub(size(G_mn),near_idx(dum));
                m = [mx,my,mz];
                % centre of the observation voxel
                r_m = ((m-1) .* d)';

                %%%[I1,I2,I3,I4] = surface_surface_kernels(W,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx);
                %%%
                %%% multiply the values by the constant coefficients
                %%% and return the sum over the 36 faces
                %%G_mn(mx,my,mz,:) = surface_surface_coeff(I1,I2,I3,I4,dx,ko_grfn); 

                % uncomment the following due to selection of the calculation
                % method

                % 1) calculate near interactions with singularity subtraction
                % method. VV for the smoothed kernel, SS for 1/R, SS for R/2
                I_V_smooth = volume_volume(W6D,X,Y,Z,Xp,Yp,Zp,ko_grfn,r_m,r_n,dx,4);

                I_S2 = surface_surface_kernels(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,2);
                %[T,I_S2] = evalc('surface_surface_kernels(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,2);');

                I_S3 = surface_surface_kernels(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,3);

                G_mn_tmp(dum,:) = I_V_smooth + I_S2 + I_S3;

                % 2) calculate near interactions with VV to SS formulas for G
                %G_mn(mx,my,mz,:) = surface_surface_kernels(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,1);
    end
    % and now properly insert the near G_mn elements in the 'G_mn' tensor. 
    for dum = 1:size(near_idx,1)
                [mx,my,mz] = ind2sub(size(G_mn),near_idx(dum));
                G_mn(mx,my,mz,:) = G_mn_tmp(dum,:);
    end
    
    disp(['  Time for computing near interactions using parallelization ::: ',num2str(toc)])

% old code, loosely parallelized    
else
    parfor mx = 1:Ls
        for my = 1:Ms
            for mz = 1:Ns

                m = [mx,my,mz];
                % centre of the observation voxel
                r_m = ((m-1) .* d)';


                % 1) calculate near interactions with singularity subtraction
            % method. volume-volume for the smoothed kernel, surface-surface for 1/R, surface-surface for R

            % non-singular, smooth kernel:  KER = 1/4/pi * ( exp(-1j*k0*R) - 1.0 ) ./R + k0^2/8/pi * R;
                I_V_smooth = volume_volume(W6D,X,Y,Z,Xp,Yp,Zp,ko_grfn,r_m,r_n,dx,4);

            % singular kernel: KER =  1/(4*pi*R) with surface-surface formulas
                I_S2 = surface_surface_kernels(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,2);

            % singular kernel: KER = -ko^2/(8*pi)*R with surface-surface formulas
                I_S3 = surface_surface_kernels(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,3);

            % take the sum of them according to eq. 17 of the report
                G_mn(mx,my,mz,:) = I_V_smooth + I_S2 + I_S3;

                % 2) calculate near interactions with VV to SS formulas for G
                %G_mn(mx,my,mz,:) = surface_surface_kernels(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,1);


           end
        end
    end
    
    disp(['  Time for computing near interactions the old way ::: ',num2str(toc)])
end



G_mn = G_mn /dx^4 * ko^2;

% Normalizing linear basis functions
% Attention::: Here has been added!!!
G_mn(:,:,:,2:3) = G_mn(:,:,:,2:3) *(1/dx);
G_mn(:,:,:,4) = G_mn(:,:,:,4) *(1/(dx^2));
G_mn(:,:,:,5) = G_mn(:,:,:,5) *(1/(dx));
G_mn(:,:,:,6:7) = G_mn(:,:,:,6:7) *(1/(dx^2));

% Self terms for sparse preconditioner

st_sparse_precon=[squeeze(G_mn(1,1,1,1)) squeeze(G_mn(1,1,1,4)) squeeze(G_mn(1,1,1,7))];

% Compute circulant tensors
tic
Gp_mn = zeros(2*L,2*M,2*N,7);
infofN = whos('Gp_mn'); 
memestimated = 2*infofN.bytes/(1024*1024); % times 2 for real2cmplx
disp(['  Memory for storing circulant tensor (MB) ::: ' , num2str(memestimated)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cube 'L'
[Gp_coeff_L] = gperiodic_coeff_nop_linJVIE('L');
% Cube 'M'
[Gp_coeff_M] = gperiodic_coeff_nop_linJVIE('M');
% Cube 'N'
[Gp_coeff_N] = gperiodic_coeff_nop_linJVIE('N');
% Cube 'LM'
[Gp_coeff_LM] = gperiodic_coeff_nop_linJVIE('LM');
% Cube 'LN'
[Gp_coeff_LN] = gperiodic_coeff_nop_linJVIE('LN');
% Cube 'MN'
[Gp_coeff_MN] = gperiodic_coeff_nop_linJVIE('MN');
% Cube 'LMN'
[Gp_coeff_LMN] = gperiodic_coeff_nop_linJVIE('LMN');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gp_mn(1:L,1:M,1:N,:) = G_mn;
%
for ii = 1:7
    % Cube 'L'
    Gp_mn(L+2:2*L,1:M,1:N,ii)         = G_mn(L:-1:2,1:M,1:N,ii) * Gp_coeff_L(ii,1);
    % Cube 'M'
    Gp_mn(1:L,M+2:2*M,1:N,ii)         = G_mn(1:L,M:-1:2,1:N,ii) * Gp_coeff_M(ii,1);
    % Cube 'N'
    Gp_mn(1:L,1:M,N+2:2*N,ii)         = G_mn(1:L,1:M,N:-1:2,ii) * Gp_coeff_N(ii,1);
    % Cube 'LM'
    Gp_mn(L+2:2*L,M+2:2*M,1:N,ii)     = G_mn(L:-1:2,M:-1:2,1:N,ii) * Gp_coeff_LM(ii,1);
    % Cube 'LN'
    Gp_mn(L+2:2*L,1:M,N+2:2*N,ii)     = G_mn(L:-1:2,1:M,N:-1:2,ii) * Gp_coeff_LN(ii,1);
    % Cube 'MN'
    Gp_mn(1:L,M+2:2*M,N+2:2*N,ii)     = G_mn(1:L,M:-1:2,N:-1:2,ii) * Gp_coeff_MN(ii,1);
    % Cube 'LMN'
    Gp_mn(L+2:2*L,M+2:2*M,N+2:2*N,ii) = G_mn(L:-1:2,M:-1:2,N:-1:2,ii) * Gp_coeff_LMN(ii,1);
end
disp(['  Time for computing circulant tensor ::: ',num2str(toc)])
% Remove unnecessary memory consuming elements
clear G_mn

% Obtain the FFTs of circulant
if (fl_no_fft == 1)
    fN_all = Gp_mn;
    infofN = whos('fN_all'); memestimated = 2*infofN.bytes/(1024*1024);
    disp(['  Two circulant tensors will be held in memory for fast sweep!']);
    disp(['  Memory for storing two circulant tensors (MB) ::: ' , num2str(memestimated)]);
else
    tic
    fN_all = fft_operator(Gp_mn);
    disp(['  Time for FFT of circulant tensor ::: ',num2str(toc)])
end

  
