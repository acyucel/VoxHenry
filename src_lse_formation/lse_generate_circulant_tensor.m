function [fN_all,st_sparse_precon] = lse_generate_circulant_tensor(dx,ko,L,M,N,fl_no_fft)

% % Code addition for running subroutine seperately
% clc; close all; clear all
% freq = 3e0; dx = 0.1e-1; EMconstants;
% L=5; M=5; N=5
% fl_comp_fullmat = 1;

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
