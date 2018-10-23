function [Zmat_large] = lse_generate_Zmat(fN_all, L, M, N)
  
% Computation of explicit Z matrix from the circulant
% 'fN_all' must be the FFT of the circulant 
% WARNING: Z matrix can be HUGE!

disp('Zmat_large matrix will be generated! ')
tic   
 
% a) Allocate system matrix
Zmat=zeros(L*M*N,L*M*N,10);
infofN = whos('Zmat');
memestimated = 2*infofN.bytes;
fprintf('Estimated memory for Zmat temporary %.3f MB.\n' , memestimated/(1024*1024));

% b) Obtain index information from a fictitious grid
grid_intcon = ones(L,M,N);
idxS = find(abs(grid_intcon(:)) > 1e-12);
clear grid_intcon
%%
% c) Perform mat-vect multiplications to get system matrix
ind_unk = 0;
[LfN,MfN,NfN,dum_d] = size(fN_all);
for kk=1:L
%for mm=1:N
   for ll=1:M
       for mm=1:N
           %for kk=1:L
           ind_unk = ind_unk + 1;
           for slct_comp=1:10
               JIn0 = zeros(L*M*N,1);
               JIn = zeros(L, M, N);
               JOut0 = zeros(L*M*N,1);
               JOut = zeros(L, M, N);
               
               JIn0(ind_unk) = 1;
               JIn(idxS) = JIn0(:);
               
               %JIn = zeros(L, M, N); JOut = zeros(L, M, N);
               %JIn(kk,ll,mm) = 1;
               %JIn(ind_unk) = 1;
               
               % Need to distinguish in case of 2D (N = 1)
               % Remark: not supporting cases when L or M is 1
               if(N > 1)
                   fJ = fftn(JIn(:,:,:),[LfN, MfN, NfN]);
               else
                   fJ = fftn(JIn(:,:),[LfN, MfN]);
               end
               Jout1 = fN_all(:,:,:,slct_comp) .* fJ;
               
               Jout1 = ifftn(Jout1);
               
               JOut(:,:,:) = Jout1(1:L,1:M,1:N);
               
               Zmat(:,ind_unk,slct_comp) = JOut(:);
           end
       end
   end
end


% forming the large Zmat without identity

Zmat_large=zeros(5*L*M*N,5*L*M*N);

infofN = whos('Zmat_large');
memestimated = 2*infofN.bytes;
fprintf('Estimated memory for Zmat_large %.3f MB.\n' , memestimated/(1024*1024));

% The blocks in large Zmat with the last indices of fN_all
% [1 0 0 2 2; 0 1 0 3 -3; 0 0 1 0 8; 4 5 0 6 9; 4 -5 7 9 10;]

dum_zeros=zeros(L*M*N,L*M*N);

Zmat_large = [Zmat(:,:,1) dum_zeros dum_zeros Zmat(:,:,2) Zmat(:,:,2); ...
   dum_zeros Zmat(:,:,1) dum_zeros Zmat(:,:,3) -Zmat(:,:,3); ...
   dum_zeros dum_zeros Zmat(:,:,1) dum_zeros Zmat(:,:,8); ...
   Zmat(:,:,4) Zmat(:,:,5) dum_zeros Zmat(:,:,6) Zmat(:,:,9); ...
   Zmat(:,:,4) -Zmat(:,:,5) Zmat(:,:,7) Zmat(:,:,9) Zmat(:,:,10);];

clear Zmat

disp(['Time for computing Zmat_large ::: ',num2str(toc)])
    
