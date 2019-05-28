function lse_sparse_precon_prepare(dx,z_real,z_imag,idxS,st_sparse_precon,nodeid_4_grnd,nodeid_4_injectcurr,Ae)
global A_inv LL UU PP QQ RR Sch_sparse slct_decomp_sch fl_cholmod A_noninv W fPSchur precond
global fl_precon_type

if (strcmp(fl_precon_type, 'no_precond') == 1)
    disp(['Not using any preconditioner'])
    return
end
    
    
%slct_decomp_sch='ldlt_decomp'; %'chol_decomp','no_decomp','ldlt_decomp','chol_decomp'
if(exist ('OCTAVE_VERSION', 'builtin') > 0)
    % Note: 'fl_cholmod' cannot be used under Octave, not supported
    %       'ldlt_decomp' has not been ported under Octave yet
    slct_decomp_sch='lu_decomp'; %'no_decomp','lu_decomp','ldlt_decomp','chol_decomp'
else 
    slct_decomp_sch='ldlt_decomp'; %'no_decomp','lu_decomp','ldlt_decomp','chol_decomp'
end

fl_cholmod = 1; % flag for using CHOLMOD in suitesparse for fast and memory-efficient LDLT factorization and inversion
fl_volt_source = 2; % symmetric voltage source implementation
fl_profile = 0; % cpu and memory profiling


fl_accu_check = 0; % flag for accuracy checking for inversion with the decompositions
tic_Assembly = tic;

% constants
num_node = size(Ae,1);
num_curr = size(Ae,2);
num_curr_one5th = num_curr/5;

% 1) Compute A_inv matrix
tic
% get the actual values of 'OneoverSigma_e'. This allows different conductivites for each voxel.
z_real_nonemptyvoxel = z_real(idxS);
% if not superconductor
if isempty(z_imag)
    diag_pulse=1./(z_real_nonemptyvoxel/dx + st_sparse_precon(1));
    diag_2Dlinear=1./(z_real_nonemptyvoxel/(6*dx) + st_sparse_precon(2));
    diag_3Dlinear=1./(z_real_nonemptyvoxel/(2*dx) + st_sparse_precon(3));
else
    z_imag_nonemptyvoxel = z_imag(idxS);
    diag_pulse=1./((z_real_nonemptyvoxel+z_imag_nonemptyvoxel)/dx + st_sparse_precon(1));
    diag_2Dlinear=1./((z_real_nonemptyvoxel+z_imag_nonemptyvoxel)/(6*dx) + st_sparse_precon(2));
    diag_3Dlinear=1./((z_real_nonemptyvoxel+z_imag_nonemptyvoxel)/(2*dx) + st_sparse_precon(3));
end
  
inds=zeros(num_curr,3);
% rows for sparse 'A_inv' formation
inds(1:num_curr,1)=[1:1:num_curr];
% columns for sparse 'A_inv' formation
inds(1:num_curr,2)=inds(1:num_curr,1);
% values for sparse 'A_inv' formation, split in three sets:
% first set is from 1 to 3/5 of num_curr, i.e. Ix, Iy, Iz
inds(1:num_curr_one5th,3)=abs(diag_pulse);
inds(num_curr_one5th+1:2*num_curr_one5th,3)=abs(diag_pulse);
inds(2*num_curr_one5th+1:3*num_curr_one5th,3)=abs(diag_pulse);
% second set is from 3/5 to 4/5 of num_curr, i.e. I2D
inds(3*num_curr_one5th+1:4*num_curr_one5th,3)=abs(diag_2Dlinear);
% third set is from 4/5 to 5/5 of num_curr, i.e. I3D
inds(4*num_curr_one5th+1:num_curr,3)=abs(diag_3Dlinear);
% now create the sparse 'A_inv'
A_inv=sparse(inds(:,1),inds(:,2),inds(:,3));
if (fl_profile == 1); disp(['Time to compute A_inv block ::: ',num2str(toc)]); end

infomem1 = whos('A_inv');
memestimated = (infomem1.bytes)/(1024*1024);
if (fl_profile == 1); disp(['Memory for A_inv block (MB)::' , num2str(memestimated)]); end

if (strcmp(fl_precon_type, 'schur_approx') == 1)
  
    % need the non-inverted A matrix (diagonal matrix approximating Z)
    A_noninv = inv(A_inv);
    
    % calculate W
    %
    % create the weights
    w = sum(abs(Ae),2);
    % and place them on the diagonal of a sparse matrix, 'W'
    W = spdiags(1./w,0,size(w,1),size(w,1));

elseif (strcmp(fl_precon_type, 'schur_gmres') == 1)
    
    Sch_comp=Ae*A_inv*Ae';
    precond = spdiags(diag(Sch_comp),0,size(Sch_comp,1),size(Sch_comp,1));

    fPSchur = @(schur_in)lse_sparse_schur_compl_multiply(schur_in, Ae, A_inv);

elseif (strcmp(fl_precon_type, 'schur_invert') == 1)
  
    % 2) Obtain Schur complement
    tic
    
    % as the rows of Ae corresponding to ground or excitation nodes have been removed,
    % obtaining what in the VoxHenry paper is called 'Ar', we can proceed with the standard
    % Schur complement from the system (but in the code 'Ar' is actually 'Ae' without the zeroed rows)
    % [A  B] = [Z -Ar']
    % [C  D]   [Ar  0 ]

    % symm. voltage source and current source will use this
    Sch_comp=Ae*A_inv*Ae';

    if (fl_profile == 1); disp(['Time to obtain Schur complement ::: ',num2str(toc)]); end

    % 3) Factorization of Schur complement with different schemes

    %ishermitian(Sch_comp)
    if (issymmetric(Sch_comp) == 0)
         if (strcmp(slct_decomp_sch, 'ldlt_decomp') == 1)
             disp('Schur complement matrix is not symmetric')
             disp('Decomposition is changed to LU')
             slct_decomp_sch = 'lu_decomp';
         elseif (strcmp(slct_decomp_sch, 'chol_decomp') == 1 )
             disp('Schur complement matrix is not symmetric')
             disp('Decomposition is changed to LU')
             slct_decomp_sch = 'lu_decomp';
         end
    end

    switch slct_decomp_sch
        
        case 'lu_decomp'
            
            tic
            [LL,UU,PP,QQ,RR] = lu(Sch_comp);
            disp(['Time to (LU) factorize Schur complement ::: ',num2str(toc)])
            
            infomem1 = whos('LL');infomem2 = whos('UU'); infomem3 = whos('PP');
            infomem4 = whos('QQ');infomem5 = whos('RR');
            memestimated = (infomem1.bytes+infomem2.bytes+infomem3.bytes+infomem4.bytes+infomem5.bytes)/(1024*1024);
            disp(['Memory for LU fact. matrices (MB)::' , num2str(memestimated)]);
            % Test for CPU time of inversion
            bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
            tic
            x_dum_lu = QQ * (UU \ (LL \ (PP * (RR \ bb))));
            disp(['Time for one inversion with LU ::: ',num2str(toc)])
            tic
            
            if (fl_accu_check == 1)
                tic
                x_dum_back = Sch_comp\bb;
                disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                csd=num2str(max(abs(x_dum_lu-x_dum_back)./abs(x_dum_back)));
                disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
            end
            
            
        case 'ldlt_decomp'
            
            if (fl_cholmod == 1)
                % factorization and inversion with CHOLMOD in Suitesparse
                tic;
                [LL,p,QQ] = ldlchol (Sch_comp);
                disp(['Time to (LDLT) factorize Schur complement - Cholmod ::: ',num2str(toc)])
                
                infomem1 = whos('LL');infomem2 = whos('p');infomem3 = whos('QQ');
                memestimated = (infomem1.bytes+infomem2.bytes+infomem3.bytes)/(1024*1024);
                disp(['Memory for LDLT fact. matrices - Cholmod (MB)::' , num2str(memestimated)]);
                
                if (fl_accu_check == 1)
                    tic
                    x_dum_back = Sch_comp\bb;
                    disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                    csd=num2str(max(abs(x_dum_ldlt-x_dum_back)./abs(x_dum_back)));
                    disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
                end
                
            else
                
                tic
                [LL,QQ,PP,RR] = ldl(Sch_comp);% option 2 PP'*RR*A*RR*PP = LL*QQ*LL'
                disp(['Time to (LDLT) factorize Schur complement ::: ',num2str(toc)])
                
                infomem1 = whos('LL');infomem2 = whos('QQ');infomem3 = whos('PP'); infomem4 = whos('RR');
                memestimated = (infomem1.bytes+infomem2.bytes+infomem3.bytes+infomem4.bytes)/(1024*1024);
                disp(['Memory for LDLT fact. matrices (MB)::' , num2str(memestimated)]);
                % Test for CPU time of inversion
                bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
                tic
                x_dum_ldlt = RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * bb)))); % for option 2
                disp(['Time for one inversion with LDLT ::: ',num2str(toc)])
                
                if (fl_accu_check == 1)
                    tic
                    x_dum_back = Sch_comp\bb;
                    disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                    csd=num2str(max(abs(x_dum_ldlt-x_dum_back)./abs(x_dum_back)));
                    disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
                end
            end


        case 'chol_decomp'
            
            tic
            [LL,g,PP] = chol(Sch_comp,'lower'); % option 3: LL'=S'AS
            if (g ~= 0)
                error ('MATLAB:posdef', 'Matrix must be positive definite.') ;
            end
            disp(['Time to (Cholesky) factorize Schur complement ::: ',num2str(toc)])
                    
            infomem1 = whos('LL');infomem2 = whos('PP');
            memestimated = (infomem1.bytes+infomem2.bytes)/(1024*1024);
            disp(['Memory for Chol. fact. matrices (MB)::' , num2str(memestimated)]);
            % Test for CPU time of inversion
            %bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
            %tic
            %x_dum_chol = PP * (LL' \ (LL \ (PP' * bb))); % for option 3
            %disp(['Time for one inversion with Chol ::: ',num2str(toc)])
            
            if (fl_accu_check == 1)
                tic
                x_dum_back = Sch_comp\bb;
                disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                csd=num2str(max(abs(x_dum_chol-x_dum_back)./abs(x_dum_back)));
                disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
            end
            
            % For large scae test, the vector permutation implementation below
            % didn't reduce the time. However, I include this here for future
            % reference

            % % Cholesky with vector permutation
            % tic
            % [LL,g,s] = chol(Sch_comp,'lower','vector');
            % disp(['Time to (Cholesky) factorize Schur complement - vector ::: ',num2str(toc)])
            % 
            % x_dum_chol3=zeros(size(x_dum_chol,1),size(x_dum_chol,2));
            % 
            % tic
            % x_dum_chol3(s) = (LL' \ (LL \ (bb(s)))); % for option 3
            % disp(['Time for one inversion with Chol - vector ::: ',num2str(toc)])
            % tic
            % 
            % rel_diff2=max(abs(x_dum_chol -x_dum_chol3)./abs(x_dum_chol))
            
            
        case 'no_decomp'
            
            infomem1 = whos('Sch_comp');
            memestimated = (infomem1.bytes)/(1024*1024);
            disp(['Memory for Schur complement matrix (MB)::' , num2str(memestimated)]);
            
            Sch_sparse = Sch_comp;
            
        otherwise
            
            error('Invalid decomposition selection for Schur complement')
            
    end
    
elseif (strcmp(fl_precon_type, 'schur_invert_original') == 1)
  
    % 2) Obtain Schur complement
    tic
    
    % as the rows of Ae corresponding to ground or excitation nodes have been zeroed,
    % the Ae*A_inv*Ae' matrix is rank deficient (singlular). DD adds some dummy equations
    % placing 1 on the diagonal of the block D (block that normally should be zero)
    % where there are zeroed rows
    % [A  B] = [Z -Ae']
    % [C  D]   [Ae  0 ]
    % 
    if (fl_volt_source == 1)
        DD = spalloc(num_node,num_node,length(nodeid_4_grnd)+length(nodeid_4_injectcurr));
    else
        DD = spalloc(num_node,num_node,length(nodeid_4_grnd));
    end

    for kk=1:length(nodeid_4_grnd)
        DD(nodeid_4_grnd(kk),nodeid_4_grnd(kk))=1;
    end

    if (fl_volt_source == 1 || fl_volt_source == 2)
        for kk=1:length(nodeid_4_injectcurr)
            DD(nodeid_4_injectcurr(kk),nodeid_4_injectcurr(kk))=1;
        end
    end

    if (fl_volt_source == 1)
        Ae_tmp=Ae; % Attention::: temporarily doubling the memory for Ae!!!
        Ae_tmp(nodeid_4_grnd,:)=0;
        Ae_tmp(nodeid_4_injectcurr,:)=0;
        Sch_comp=DD+Ae_tmp*A_inv*Ae';
        clear Ae_tmp
    else % symm. voltage source and current source will use this
        Sch_comp=DD+Ae*A_inv*Ae';
        %Sch_comp=Ae*A_inv*Ae';
    end
    if (fl_profile == 1); disp(['Time to obtain Schur complement ::: ',num2str(toc)]); end


    % 3) Factorization of Schur complement with different schemes

    %ishermitian(Sch_comp)
     if (issymmetric(Sch_comp) == 0)
         if (strcmp(slct_decomp_sch, 'ldlt_decomp') == 1)
             disp('Schur complement matrix is not symmetric')
             disp('Decomposition is changed to LU')
             slct_decomp_sch = 'lu_decomp';
         elseif (strcmp(slct_decomp_sch, 'chol_decomp') == 1 )
             disp('Schur complement matrix is not symmetric')
             disp('Decomposition is changed to LU')
             slct_decomp_sch = 'lu_decomp';
         end
     end

    switch slct_decomp_sch
        
        case 'lu_decomp'
            
            tic
            [LL,UU,PP,QQ,RR] = lu(Sch_comp);
            disp(['Time to (LU) factorize Schur complement ::: ',num2str(toc)])
            
            infomem1 = whos('LL');infomem2 = whos('UU'); infomem3 = whos('PP');
            infomem4 = whos('QQ');infomem5 = whos('RR');
            memestimated = (infomem1.bytes+infomem2.bytes+infomem3.bytes+infomem4.bytes+infomem5.bytes)/(1024*1024);
            disp(['Memory for LU fact. matrices (MB)::' , num2str(memestimated)]);
            % Test for CPU time of inversion
            bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
            tic
            x_dum_lu = QQ * (UU \ (LL \ (PP * (RR \ bb))));
            disp(['Time for one inversion with LU ::: ',num2str(toc)])
            tic
            
            if (fl_accu_check == 1)
                tic
                x_dum_back = Sch_comp\bb;
                disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                csd=num2str(max(abs(x_dum_lu-x_dum_back)./abs(x_dum_back)));
                disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
            end
            
            
        case 'ldlt_decomp'
            
            if (fl_volt_source == 1)
                error('Sch comp not symmetric for voltage source - LDLT decomp does not work - change decomposition')
            end
            
            if (fl_cholmod == 1)
                % factorization and inversion with CHOLMOD in Suitesparse
                tic;
                [LL,p,QQ] = ldlchol (Sch_comp);
                disp(['Time to (LDLT) factorize Schur complement - Cholmod ::: ',num2str(toc)])
                
                infomem1 = whos('LL');infomem2 = whos('p');infomem3 = whos('QQ');
                memestimated = (infomem1.bytes+infomem2.bytes+infomem3.bytes)/(1024*1024);
                disp(['Memory for LDLT fact. matrices - Cholmod (MB)::' , num2str(memestimated)]);
                
                %x_dum_ldlt=zeros(size(Sch_comp,1),1);
                %bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
                
                %tic;
                %x_dum_ldlt(QQ) = ldlsolve (LL,bb(QQ));
                %disp(['Time for one inversion with LDLT - Cholmod ::: ',num2str(toc)])
                
                if (fl_accu_check == 1)
                    tic
                    x_dum_back = Sch_comp\bb;
                    disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                    csd=num2str(max(abs(x_dum_ldlt-x_dum_back)./abs(x_dum_back)));
                    disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
                end
                
            else
                
                tic
                [LL,QQ,PP,RR] = ldl(Sch_comp);% option 2 PP'*RR*A*RR*PP = LL*QQ*LL'
                disp(['Time to (LDLT) factorize Schur complement ::: ',num2str(toc)])
                
                infomem1 = whos('LL');infomem2 = whos('QQ');infomem3 = whos('PP'); infomem4 = whos('RR');
                memestimated = (infomem1.bytes+infomem2.bytes+infomem3.bytes+infomem4.bytes)/(1024*1024);
                disp(['Memory for LDLT fact. matrices (MB)::' , num2str(memestimated)]);
                % Test for CPU time of inversion
                bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
                tic
                x_dum_ldlt = RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * bb)))); % for option 2
                disp(['Time for one inversion with LDLT ::: ',num2str(toc)])
                
                if (fl_accu_check == 1)
                    tic
                    x_dum_back = Sch_comp\bb;
                    disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                    csd=num2str(max(abs(x_dum_ldlt-x_dum_back)./abs(x_dum_back)));
                    disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
                end
            end


        case 'chol_decomp'
            
            if (fl_volt_source == 1)
                error('Sch comp not symmetric for voltage source - Cholesky decomp does not work - change decomposition')
            end
            
            tic
            [LL,g,PP] = chol(Sch_comp,'lower'); % option 3: LL'=S'AS
            if (g ~= 0)
                error ('MATLAB:posdef', 'Matrix must be positive definite.') ;
            end
            disp(['Time to (Cholesky) factorize Schur complement ::: ',num2str(toc)])
                    
            infomem1 = whos('LL');infomem2 = whos('PP');
            memestimated = (infomem1.bytes+infomem2.bytes)/(1024*1024);
            disp(['Memory for Chol. fact. matrices (MB)::' , num2str(memestimated)]);
            % Test for CPU time of inversion
            %bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
            %tic
            %x_dum_chol = PP * (LL' \ (LL \ (PP' * bb))); % for option 3
            %disp(['Time for one inversion with Chol ::: ',num2str(toc)])
            
            if (fl_accu_check == 1)
                tic
                x_dum_back = Sch_comp\bb;
                disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                csd=num2str(max(abs(x_dum_chol-x_dum_back)./abs(x_dum_back)));
                disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
            end
            
            % For large scae test, the vector permutation implementation below
            % didn't reduce the time. However, I include this here for future
            % reference

            % % Cholesky with vector permutation
            % tic
            % [LL,g,s] = chol(Sch_comp,'lower','vector');
            % disp(['Time to (Cholesky) factorize Schur complement - vector ::: ',num2str(toc)])
            % 
            % x_dum_chol3=zeros(size(x_dum_chol,1),size(x_dum_chol,2));
            % 
            % tic
            % x_dum_chol3(s) = (LL' \ (LL \ (bb(s)))); % for option 3
            % disp(['Time for one inversion with Chol - vector ::: ',num2str(toc)])
            % tic
            % 
            % rel_diff2=max(abs(x_dum_chol -x_dum_chol3)./abs(x_dum_chol))
            
            
        case 'no_decomp'
            
            infomem1 = whos('Sch_comp');
            memestimated = (infomem1.bytes)/(1024*1024);
            disp(['Memory for Schur complement matrix (MB)::' , num2str(memestimated)]);
            
            Sch_sparse = Sch_comp;
            
        otherwise
            
            error('Invalid decomposition selection for Schur complement')
            
    end
end

Time_Assembly = toc(tic_Assembly);
disp(['Time to prepare sparse precon ::: ',num2str(Time_Assembly)])
