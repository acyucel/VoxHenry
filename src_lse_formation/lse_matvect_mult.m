function [JOut_full]=lse_matvect_mult(JIn0, fN_all, Ae, OneoverMc, dx, freq, idx, nodeid_4_grnd,nodeid_4_injectcurr)
global fl_precon_type

% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

fl_volt_source = 2; %symmetric voltage source
fl_gpu = 0;
fl_profile = 0;

tic
% constants
mu = 4*pi*1e-7;
co = 299792458; 
eo = 1/co^2/mu;
omega = 2*pi*freq;
oneoverjomegaeo=1/(j*omega*eo);

num_node=size(Ae,1);
num_curr=size(Ae,2);

% fft dimensions
[LfN, MfN, NfN, ~] = size(fN_all);

% domain dimensions
[L, M, N] = size(OneoverMc);

%GPU_flag = 0;

if (fl_gpu == 1)
    % allocate space
    JIn = gpuArray.zeros(L, M, N, 5);    
    JOut = gpuArray.zeros(L, M, N, 5);
    % send to gpu, and translate from local (idx) to global (L,M,N) coordinates
    JIn(idx) = gpuArray(JIn0(1:num_curr));
else

    % allocate space
    JIn = zeros(L, M, N, 5);
    JOut = zeros(L, M, N, 5);
   
    % translate from local (idx) to global (L,M,N) coordinates
    JIn(idx) = JIn0(1:num_curr);
end

JOut_full = zeros(num_curr+num_node,1);
JOut_full_in = zeros(num_curr+num_node,1);

%%% allocate space
%JIn = zeros(L, M, N, 3);
%JOut = zeros(L, M, N, 3);
%JOut_full = zeros(num_curr+num_node,1);

% translate from local (idx) to global (L,M,N) coordinates
%JIn(idx) = JIn0(1:num_curr);

% ---------------------------------------------------------------------
% apply fft and mv-op for each of the components of JIn
% ---------------------------------------------------------------------

% x component of JIn, store contribution on 3 components of Jout
fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
Jout1 = fN_all(:,:,:,1) .* fJ; % Gxx*Jx
Jout4 = -fN_all(:,:,:,2) .* fJ; % G2dx*Jx
Jout5 = Jout4; % G3dx*Jx

% y component of JIn, add contribution on 3 components of Jout
fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
Jout2 = fN_all(:,:,:,1) .* fJ; % Gyy*Jy
Jout4 = Jout4 - fN_all(:,:,:,3) .* fJ; % G2dy*Jy - (-(y-yc))
Jout5 = Jout5 + fN_all(:,:,:,3) .* fJ; % G3dy*Jy - (y-yc)

% z component of JIn, store contribution on 2 components of Jout
fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
Jout3 = fN_all(:,:,:,1) .* fJ; % Gzz*Jz
Jout5 = Jout5 - fN_all(:,:,:,5) .* fJ; % G3dz*Jz

% 2d component of JIn, add contribution on 4 components of Jout
fJ = fftn(JIn(:,:,:,4),[LfN, MfN, NfN]);
Jout1 = Jout1 + fN_all(:,:,:,2) .* fJ; % Gx2d*J2d
Jout2 = Jout2 + fN_all(:,:,:,3) .* fJ; % Gy2d*J2d - (-(y-yc))
Jout4 = Jout4 + fN_all(:,:,:,4) .* fJ; % G2d2d*J2d
Jout5 = Jout5 + fN_all(:,:,:,6) .* fJ; % G3d2d*J2d

% 3d component of JIn, add contribution on 3 components of Jout
fJ = fftn(JIn(:,:,:,5),[LfN, MfN, NfN]);
Jout1 = Jout1 + fN_all(:,:,:,2) .* fJ; % Gx3d*J3d
Jout2 = Jout2 - fN_all(:,:,:,3) .* fJ; % Gy3d*J3d - (y-yc)
Jout3 = Jout3 + fN_all(:,:,:,5) .* fJ; % Gz3d*J3d
Jout4 = Jout4 + fN_all(:,:,:,6) .* fJ; % G2d3d*J3d
Jout5 = Jout5 + fN_all(:,:,:,7) .* fJ; % G3d3d*J3d



% apply ifft and add identity term
Jout1 = ifftn(Jout1);
JOut(:,:,:,1) = (1/dx) .* OneoverMc .* JIn(:,:,:,1) - Jout1(1:L,1:M,1:N);
Jout2 = ifftn(Jout2);
JOut(:,:,:,2) = (1/dx) .* OneoverMc .* JIn(:,:,:,2) - Jout2(1:L,1:M,1:N);
Jout3 = ifftn(Jout3);
JOut(:,:,:,3) = (1/dx) .* OneoverMc .* JIn(:,:,:,3) - Jout3(1:L,1:M,1:N);
Jout4 = ifftn(Jout4);
JOut(:,:,:,4) = (dx/6/(dx^2)) .* OneoverMc .* JIn(:,:,:,4) - Jout4(1:L,1:M,1:N);
Jout5 = ifftn(Jout5);
JOut(:,:,:,5) = (dx/2/(dx^2)) .* OneoverMc .* JIn(:,:,:,5) - Jout5(1:L,1:M,1:N);

% multiply by 1/(jweps0)

JOut(:,:,:,:) = JOut(:,:,:,:) * oneoverjomegaeo; 

% -------------------------------------------------------------------------
% Return local coordinates related to material positions
% -------------------------------------------------------------------------

if (fl_gpu == 1)
    % get from GPU
    JOut = gather(JOut(idx));
    % clear gpu data
    clear JIn; clear Jout1; clear Jout2; clear Jout3; clear fJ;
else
    JOut = JOut(idx);
end

%JOut = JOut(idx);

JOut_full(1:num_curr) = JOut; 

if(fl_profile == 1); disp(['Time for matvect - fft part::: ',num2str(toc)]); end;

% ---------------------------------------------------------------------
% Adding contributions due to nodal incidence matrix
% ---------------------------------------------------------------------
tic
% Perform multiplications without assigning to dum_block

JOut_full(1:num_curr) = JOut_full(1:num_curr) - (Ae'*JIn0(num_curr+1:num_curr+num_node)) ;

JOut_full(num_curr+1:num_curr+num_node) = Ae*JIn0(1:num_curr);

% this is needed only in case of the original code Schur inversion method, because in this
% case instead of removing the empty rows of Ae, they have been zeroed, so there is the need
% to add dummy equations to the system (see comments in 'lse_sparse_precond_prepare.m'
% relevant to 'DD' matrix)
if ( strcmp(fl_precon_type, 'schur_invert_original') == 1 )
    % For "well-conditioning"
    JOut_full(num_curr+nodeid_4_grnd) = JIn0(num_curr+nodeid_4_grnd);

    if(fl_volt_source == 1 || fl_volt_source == 2)
        JOut_full(num_curr+nodeid_4_injectcurr) = JIn0(num_curr+nodeid_4_injectcurr);
    end
end

if(fl_profile == 1); disp(['Time for matvect - Ae matrices part::: ',num2str(toc)]); end

% ---------------------------------------------------------------------
% Sparse preconditioner [E F; G H]
% ---------------------------------------------------------------------
%tic
%if ( (strcmp(fl_precon_type, 'no_precond') == 0) & (strcmp(fl_precon_type, 'test_no_precond') == 0) )
%    [JOut_full]=lse_sparse_precon_multiply(JOut_full,Ae,nodeid_4_grnd,nodeid_4_injectcurr);
%end
%if(fl_profile == 1); disp(['Time for matvect - sparse preconditioner part::: ',num2str(toc)]); end

% JOut_full_in = JOut_full;
% % block E contribution
% JOut_full(1:num_curr) = A_inv2*JOut_full_in(1:num_curr)+A_inv2*(-Ae')*...
%     QQ * (UU \ (LL \ (PP * (RR \ (Ae*A_inv2*JOut_full_in(1:num_curr))))));
% 
% % block F contribution
% JOut_full(1:num_curr) = JOut_full(1:num_curr)...
%     + A_inv2 * (Ae') * QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));
% 
% 
% % block G contribution
% JOut_full(num_curr+1:num_curr+num_node) = ...
%     -QQ * (UU \ (LL \ (PP * (RR \ (Ae*A_inv2*JOut_full_in(1:num_curr))))));
% % block H contribution
% JOut_full(num_curr+1:num_curr+num_node) = JOut_full(num_curr+1:num_curr+num_node)...
%     +QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));

% typing one dot for each iteratative solution multiplication step
fprintf ('.') ;

