function [fN_all,st_sparse_precon] = lse_generate_circulant_tensor(dx,ko,L,M,N,fl_no_fft)

% Code addition for running subroutine seperately
% clc; close all; clear all
% freq = 3e0; dx = 0.1e-1; EMconstants;
% L=3; M=3; N=3
% fl_comp_fullmat = 1;

% if 'lse_enhanced_parallel' == 1, unroll the for loops to improve
% parallelization for near interactions calculation
lse_enhanced_parallel = 1;
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
Np_1D_far_V    = 3;
% number of 1D points for medium volume-volume integrals
Np_1D_medium_V = 5;
% number of 1D points for near surface-surface and singular integrals 
Np_1D_near_S   = 12;

% Set region of medium distances cells for higher order quadrature
n_medium = 5; 

% centre of the source voxel - reference cell
r_n = [0 0 0]';
% distance vector
dy = dx; dz = dx; 
d = [dx,dy,dz];

% 2b) Compute interaction between far distance cells - the volume-volume integrals
G_mn = zeros(L,M,N,10);
infofN = whos('G_mn'); memestimated = 2*infofN.bytes/(1024*1024);
disp(['Memory for temporarily storing G_mn (MB) ::: ' , num2str(memestimated)]);

Np = Np_1D_far_V;
% 6D Clenshaw - Curtis weights and points
[W6D,X,Y,Z,Xp,Yp,Zp] = weights_points(Np,6);
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
disp(['Time for computing far interactions ::: ',num2str(toc)])

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
if(lse_exploit_symmetry == 0 || n_mediumL != n_mediumM || n_mediumL != n_mediumN ) 
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
    disp(['Time for computing medium interactions ::: ',num2str(toc)])
else
    % test to compute medium interactions exploiting symmetry
    % remark: volume must be cubic, i.e. 'mediumL' = 'mediumM' = 'mediumN'
    tic
    G_mn_test = zeros(n_mediumL,n_mediumL,n_mediumL,10);
    for mx = 1:n_mediumL
        for my = 1:n_mediumM
            % could further exploit symmetry
            for mz = 1:n_mediumN
            %for mz = 1:my
                
                m = [mx,my,mz];
                % centre of the observation voxel
                r_m = ((m-1) .* d)';
                
                K = volume_volume_sym(W6D,X,Y,Z,Xp,Yp,Zp,ko_grfn,r_m,r_n,dx,1,1);
                
                % this is only half (1/8 of complete angle around the x-axis);
                % must copy to the other half, to cover the full quadrant
                
                % Gx,x
                G_mn_test(mx,my,mz,1) = K(1);
                % G2D,x
                G_mn_test(mx,my,mz,2) = -K(2);
                % G2D,y
                G_mn_test(mz,mx,my,3) = K(2);
                % Gx,2D
                G_mn_test(mx,my,mz,4) = K(2);
                % Gy,2D
                G_mn_test(mz,mx,my,5) = -K(2);
                % G2D,2D
                % linear(x)-linear(x') part
                G_mn_test(mx,my,mz,6) = G_mn_test(mx,my,mz,6) + K(3);
                % linear(y)-linear(y') part
                G_mn_test(mz,mx,my,6) = G_mn_test(mz,mx,my,6) + K(3);
                % Gz,3D
                G_mn_test(my,mz,mx,7) = -2*K(2);
                % G3D,z
                G_mn_test(my,mz,mx,8) = 2*K(2);
                % G2D,3D
                % linear(x)-linear(x') part
                G_mn_test(mx,my,mz,9) = G_mn_test(mx,my,mz,9) + K(3);
                % linear(y)-linear(y') part
                G_mn_test(mz,mx,my,9) = G_mn_test(mz,mx,my,9) - K(3);
                % G3D,3D
                % linear(x)-linear(x') part
                G_mn_test(mx,my,mz,10) = G_mn_test(mx,my,mz,10) + K(3);
                % linear(y)-linear(y') part
                G_mn_test(mz,mx,my,10) = G_mn_test(mz,mx,my,10) + K(3);
                % linear(z)-linear(z') part
                G_mn_test(my,mz,mx,10) = G_mn_test(my,mz,mx,10) + 4*K(3);
                
            end
        end
    end
    disp(['Time for computing medium interactions exploiting symmetry ::: ',num2str(toc)])
    % delete following line and use G_mn directly when fully tested
    G_mn(1:n_mediumL,1:n_mediumL,1:n_mediumL,:) = G_mn_test(1:n_mediumL,1:n_mediumL,1:n_mediumL,:);
end

% 2e) Compute interaction between near distance cells - the surface-surface integrals

Np = Np_1D_near_S;
% 4D Clenshaw - Curtis points and weights
[W4D,A,B,C,D] = weights_points(Np,4);

Ls = 2;Ms = 2;Ns = 2;
% adjust domain for singular integrals calculations
if L <2; Ls = 1; end
if M <2; Ms = 1; end
if N <2; Ns = 1;end

tic
[I1_co, I2_co, I3_co, I4_co]  = surface_surface_coeff(dx,ko_grfn);

if (lse_enhanced_parallel == 1)
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
    G_mn_tmp = zeros(Ls * Ms * Ns, 10);
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
% old code, loosely parallelized    
else
    parfor mx = 1:Ls
        for my = 1:Ms
            for mz = 1:Ns

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

                G_mn(mx,my,mz,:) = I_V_smooth + I_S2 + I_S3;

                % 2) calculate near interactions with VV to SS formulas for G
                %G_mn(mx,my,mz,:) = surface_surface_kernels(I1_co,I2_co,I3_co,I4_co,W4D,A,B,C,D,Np,ko_grfn,r_m,r_n,m,dx,1);


           end
        end
    end
    
end

disp(['Time for computing near interactions ::: ',num2str(toc)])

G_mn = G_mn /dx^4 * ko^2;

% Normalizing linear basis functions
% Attention::: Here has been added!!!
G_mn(:,:,:,2:5) = G_mn(:,:,:,2:5) *(1/dx);
G_mn(:,:,:,6) = G_mn(:,:,:,6) *(1/(dx^2));
G_mn(:,:,:,7:8) = G_mn(:,:,:,7:8) *(1/(dx));
G_mn(:,:,:,9:10) = G_mn(:,:,:,9:10) *(1/(dx^2));

% Self terms for sparse preconditioner

st_sparse_precon=[squeeze(G_mn(1,1,1,1)) squeeze(G_mn(1,1,1,6)) squeeze(G_mn(1,1,1,10))];

% Compute circulant tensors
tic
Gp_mn = zeros(2*L,2*M,2*N ,10);
infofN = whos('Gp_mn'); memestimated = 2*infofN.bytes/(1024*1024); % times 2 for real2cmplx
disp(['Memory for storing circulant tensor (MB) ::: ' , num2str(memestimated)]);
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
for ii = 1:10
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
disp(['Time for computing circulant tensor ::: ',num2str(toc)])
% Remove unnecessary memory consuming elements
clear G_mn

% Obtain the FFTs of circulant
if (fl_no_fft == 1)
    fN_all = Gp_mn;
    infofN = whos('fN_all'); memestimated = 2*infofN.bytes/(1024*1024);
    disp(['Two circulant tensors will be held in memory for fast sweep!']);
    disp(['Memory for storing two circulant tensors (MB) ::: ' , num2str(memestimated)]);
else
    tic
    fN_all = fft_operator(Gp_mn);
    disp(['Time for FFT of circulant tensor ::: ',num2str(toc)])
end

% %% Uncomment the following if you'd like to have Z matrix
% % Compute Zmatrix from circulant and compare
% if (fl_comp_fullmat == 1)
%     
%     disp('Zmat_large matrix will be generated! ')
%     tic
%     % a) Allocate system matrix
%     Zmat=zeros(L*M*N,L*M*N,10);
%     infofN = whos('Zmat');
%     memestimated = 2*infofN.bytes;
%     fprintf('Estimated memory for Zmat temporary %.3f MB.\n' , memestimated/(1024*1024));
%     
%     % b) Obtain index information from a fictitious grid
%     grid_intcon = ones(L,M,N);
%     idxS = find(abs(grid_intcon(:)) > 1e-12);
%     clear grid_intcon
%     %%
%     % c) Perform mat-vect multiplications to get system matrix
%     ind_unk = 0;
%     [LfN,MfN,NfN,dum_d] = size(fN_all);
%     for kk=1:L
%     %for mm=1:N
%         for ll=1:M
%             for mm=1:N
%                 %for kk=1:L
%                 ind_unk = ind_unk + 1;
%                 for slct_comp=1:10
%                     JIn0 = zeros(L*M*N,1);
%                     JIn = zeros(L, M, N);
%                     JOut0 = zeros(L*M*N,1);
%                     JOut = zeros(L, M, N);
%                     
%                     JIn0(ind_unk) = 1;
%                     JIn(idxS) = JIn0(:);
%                     
%                     %JIn = zeros(L, M, N); JOut = zeros(L, M, N);
%                     %JIn(kk,ll,mm) = 1;
%                     %JIn(ind_unk) = 1;
%                     
%                     fJ = fftn(JIn(:,:,:),[LfN, MfN, NfN]);
%                     Jout1 = fN_all(:,:,:,slct_comp) .* fJ;
%                     
%                     Jout1 = ifftn(Jout1);
%                     
%                     JOut(:,:,:) = Jout1(1:L,1:M,1:N);
%                     
%                     Zmat(:,ind_unk,slct_comp) = JOut(:);
%                 end
%             end
%         end
%     end
%     
%     
%     % forming the large Zmat without identity
%     
%     Zmat_large=zeros(5*L*M*N,5*L*M*N);
%     
%     infofN = whos('Zmat_large');
%     memestimated = 2*infofN.bytes;
%     fprintf('Estimated memory for Zmat_large %.3f MB.\n' , memestimated/(1024*1024));
%     
%     % The blocks in large Zmat with the last indices of fN_all
%     % [1 0 0 2 2; 0 1 0 3 -3; 0 0 1 0 8; 4 5 0 6 9; 4 -5 7 9 10;]
%     
%     dum_zeros=zeros(L*M*N,L*M*N);
%     
%     Zmat_large = [Zmat(:,:,1) dum_zeros dum_zeros Zmat(:,:,2) Zmat(:,:,2); ...
%         dum_zeros Zmat(:,:,1) dum_zeros Zmat(:,:,3) -Zmat(:,:,3); ...
%         dum_zeros dum_zeros Zmat(:,:,1) dum_zeros Zmat(:,:,8); ...
%         Zmat(:,:,4) Zmat(:,:,5) dum_zeros Zmat(:,:,6) Zmat(:,:,9); ...
%         Zmat(:,:,4) -Zmat(:,:,5) Zmat(:,:,7) Zmat(:,:,9) Zmat(:,:,10);];
%     
%     clear Zmat
%     
%     disp(['Time for computing Zmat_large ::: ',num2str(toc)])
%     
% else
%     
%     Zmat_large = [];
%     
% end
