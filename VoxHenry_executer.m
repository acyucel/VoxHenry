clc; close all; clear all; format long e;
global fl_precon_type

% for debug only
global A_inv LL UU PP QQ RR Sch_sparse slct_decomp_sch fl_cholmod A_noninv W

disp('------------------------------------------------------------------')
disp('VoxHenry: Inductance Extraction Simulator for Voxelized Geometries')
disp('                                                                  ')
disp('          by Abdulkadir C. Yucel and Jacob K. White (MIT)         ')
disp('                                                                  ')
disp('                 Empowered by modules/ideas from                  ')
disp('   Athanasios G. Polimeridis and Ioannis P. Georgakis (Skoltech)  ')
disp('  Hakan Bagci (KAUST), Enrico Di Lorenzo (FastFieldSolvers S.R.L.)')
disp('------------------------------------------------------------------')

% -------------------------------------------------------------------------
%                  Add the Current Path to Workspace
% -------------------------------------------------------------------------

pre_define_the_path_for_folders

% -------------------------------------------------------------------------
%                  Simulation Parameters
% -------------------------------------------------------------------------

% input file
disp('Select input file:::');
inputfilelist = dir(['Input_files', filesep, '*.vhr']);
for filenum = 1:size(inputfilelist,1)
    disp([num2str(filenum), ' : ', inputfilelist(filenum).name]);
end
filenum = input('Open file number: ');
fileinname = inputfilelist(filenum).name;

er = 0;  % epsilon_r of conductors
inner_it = 100; outer_it = 10; tol=1e-12; % iterative solver inputs
prectol = 1e-1; % tolerance used in Schur GMRES inversion by 'schur_gmres' preconditioner
% use or do not use preconditioner
% valid values are: 'no_precond', 'schur_invert', 'schur_approx'
%fl_precon_type = 'no_precond';
%fl_precon_type = 'schur_approx';
fl_precon_type = 'schur_invert';
%fl_precon_type = 'schur_gmres';
%fl_precon_type = 'schur_invert_original';

% plotting
plot_currents_post_proc = 1; % if 1, plot the current densities in 3D
plot_option=1; % see the options of plotting in Visualization part
freq_curr_plot=2.5e9; % frequency for plotting currents


% -------------------------------------------------------------------------
%                  Inputs for Simulation
% -------------------------------------------------------------------------

disp(['Reading input file: ', fileinname]);
% read data from input file
% 'sigma_e' is a LxMxN array of conductivities. Zero means empty voxel, if not superconductor.
%           If superconductor, also 'lambdaL' must be zero to indicate an empty voxel.
% 'lambdaL' is a LxMxN array of London penetration depths, in case of superconductors.
%           If zero, and sigma_e is zero, means empty voxel. If there isn't any superconductor at all,
%           'lambdaL' is empty ( lambdaL = [] ) 
% 'freq' is the array of required simulation frequencies
% 'dx' is the voxel side dimension, in meters
% 'pnt_lft' is the cell array of port node relative positions, positive side 
% 'pnt_rght' is the cell array of port node relative positions, negative side 
% 'pnt_well_cond' is the cell array of grounded nodes relative positions 
[sigma_e, lambdaL, freq, dx, num_ports, pnt_lft, pnt_rght, pnt_well_cond] = pre_input_file(['Input_files', filesep, fileinname]);

% -------------------------------------------------------------------------
%                         Initialize stuff
% -------------------------------------------------------------------------

% sort and arrange frequency array
num_freq = length(freq);
if (issorted(freq) == 0) % not sorted
    freq=sort(freq)
end
freq_all = freq;
freq = freq(1); % currently do everything for the lowest freq

% generate EM constants (also frequency 'freq' dependent', e.g. 'omega')
EMconstants

% -------------------------------------------------------------------------
%                 Define EM Vars/Constants and Domain Parameters
% -------------------------------------------------------------------------

pre_define_structure_params

% output to screen a summary of the parameters
pre_print_out_inputs_generate_consts
    
% -------------------------------------------------------------------------
%                  Obtain Nodal Incidence Matrix
% -------------------------------------------------------------------------

disp('-----------------------------------------------------')
disp(['Generating panel IDs and Ae matrix...']);
tini = tic;

% Get the voxel2nodes matrix (each row is a non-empty panel, each of the six columns
% is a node ID in the order -x/+x/-y/+y/-z/+z)
%
% 'all_panels_ids' is the array of panel IDs, indexed per non-empty voxel. Voxel order
%                  is the MatLab/Octave native one for non-empty voxels of 'sigma_e'.
%                  Panel order is -x/+x/-y/+y/-z/+z. E.g. all_panels_ids(3,2) stores 
%                  the ID of the panel on the positive y face of the non-empty voxel number 3
[all_panels_ids, num_nodes] = lse_generate_panelIDs(idxS, L, M, N); 

[Ae, Ae_only_leaving, Ae_only_entering_bndry] = lse_compute_Ae_matrix(idxS, all_panels_ids, num_nodes);

tend = toc(tini);
disp(['Time for generating Ae mat & finding IDs of port nodes::: ' ,num2str(tend)]);
disp('-----------------------------------------------------')

sim_CPU_pre(1)=toc(tini); % CPU time for Ae

% ------------------------------------------------------------------------
%              Precomputation of LSE data
% ------------------------------------------------------------------------

tinisim = tic;
disp('-----------------------------------------------------')
disp(['Precomputing LSE data structures...'])

%
% circulant tensor
%

disp([' Generating circulant tensor...']);
tini = tic;

% no need to delay computation of FFT. FFT is a linear operator, and as the frequency-dependent
% part is not in the circulant tensor, we can just multiply later on by 'ko^2'
fl_no_fft=0;
% note: must still multiply 'fN_all2' and 'st_sparse_precon2' by 'ko^2' ('ko' is frequency-dependent)
% ( here we set ko = 1 in the the second parameter when calling 'lse_generate_circulant_tensor',
% so 'lse_generate_circulant_tensor' will not multiply by the actual 'ko')
[fN_all2,st_sparse_precon2] = lse_generate_circulant_tensor(dx,1,L,M,N,fl_no_fft);


tend = toc(tini);
disp([' Total time for getting circulant tensor ::: ' ,num2str(tend)]);
   
%
% port nodes and right-hand side(s)
%
 
disp([' Generating right hand side(s)...']);
tinix = tic;

% sizes of the arrays
%nodes_lft_len = sum(cellfun('size', pnt_lft, 1));
%nodes_rght_len = sum(cellfun('size', pnt_rght, 1));

% init vectors
nodeid_4_injectcurr = [];
nodeid_4_grnd = [];
num_node = size(Ae,1);
num_curr = size(Ae,2);
rhs_vect=zeros(num_node+num_curr, num_ports);

for port_no=1:num_ports
    % find the port node list for the port 'port_no'
    [port_nodeid_4_grnd, port_nodeid_4_injectcurr] = lse_assign_exc_grnd_nodes(idxS, sigma_e, all_panels_ids, port_no, pnt_lft, pnt_rght, pnt_well_cond);
    % and accumulate it
    nodeid_4_injectcurr = [nodeid_4_injectcurr; port_nodeid_4_injectcurr];
    nodeid_4_grnd = [nodeid_4_grnd; port_nodeid_4_grnd];
    
    % compute the right hand side vector, based on the relevant port node IDs
    %
    % sum the entries of excitation nodes
    exc_vect=sum(-Ae(port_nodeid_4_injectcurr, :),1);
    % just assign to the right hand side vector the first 'num_curr' elements, the rest is zero
    [rhs_vect(1:num_curr, port_no)] = full(exc_vect'); 
end

% Zeros rows and columns of Ae corresponding to ground and excitation nodes
Ae(nodeid_4_grnd,:)=0;
Ae(nodeid_4_injectcurr,:)=0;

% only the original code needs 'Ae' with rows of zeros corresponding to the ports.
% For all other cases, remove from 'Ae' and from 'rhs_vect' the port nodes rows.
if (strcmp(fl_precon_type, 'schur_invert_original') == 0)
    % length of rhs_vect is reduced of ('Ae' rows - 'Ar' rows) elements, where Ar is Ae without the rows 
    % corresponding to the ground or excitation nodes (anyway after the first size(Ae,2) elements
    % corresponding to the input port voltages, the 'rhs_vect' is all zeros)
    rhs_vect = rhs_vect(1:size(rhs_vect)-size(nodeid_4_grnd,1)-size(nodeid_4_injectcurr,1), :);
    % now remove the empty Ae rows. This should actually be 'Ar' as per VoxHenry paper
    Ae = Ae(any(Ae,2),:);
end

tendx = toc(tinix);
disp([' Total time for generating right hand side(s) ::: ' ,num2str(tendx)]);

disp(['Done... Precomputing LSE data structures'])
disp('-----------------------------------------------------')
sim_CPU_pre(2)=toc(tinisim); % CPU time for circulant+rhs


% ------------------------------------------------------------------------
%              Solving LSE
% ------------------------------------------------------------------------


disp('-----------------------------------------------------')
disp(['Solving LSEs ...'])

Y_mat=zeros(num_ports,num_ports,num_freq);
Y_mat2=zeros(num_ports,num_ports,num_freq);
Z_mat=zeros(num_ports,num_ports,num_freq);
R_jL_mat=zeros(num_ports,num_ports,num_freq);
for freq_no=1:num_freq
    tinisim = tic;

    freq = freq_all(freq_no);
    if (freq < 1e6)
        tol = 1e-12;
    else
        tol = 1e-8;
    end

    EMconstants
    disp('-----------------------------------------------------')
    disp(['Simulation for frequency : ',num2str(freq),' started! ', 'freq pnt: ',num2str(freq_no), ' / ', num2str(num_freq)])
    disp('-----------------------------------------------------')
    
    % setting new constitutive parameters for new freq
    Mr = epsilon_r - 1j*sigma_e/(eo*omega); % complex relative permittivity
    Mc = Mr - 1.0; % susceptibility
    OneoverMc = 1.0 ./ Mc; % one over susceptibility
    
    % if 'lambdaL' contains no element, then we have a normal conductor
    if isempty(lambdaL)
        # In case of normal conductors, we just need 1/sigma_e. This parameter is not
        # frequency-dependent
        z_real = 1.0 ./ sigma_e; # one over sigma
        z_imag = [];
    else      
        # General formula for superconductors:
        #
        # 1       sigma_e*(omega*mu*lambdaL^2)^2              mu*lambdaL^2
        # ----- = -------------------------------- + 1j*omega*------------------------------- = z_real + z_imag
        # sigma   (sigma_e*omega*mu*lambdaL^2)^2+1            (sigma_e*omega*mu*lambaL^2)^2+1
        %
        den = (omega*mu*sigma_e.*(lambdaL.^2)).^2 + 1;
        z_real = (sigma_e.*((omega*mu*lambdaL.^2).^2)) ./ den;
        z_imag = (1j*omega*mu*lambdaL.^2) ./ den;
        #
        # Special case: standard conductor
        # in our case, we use lambaL = 0 for a conductor that is only conducting and not superconducting.
        # As actually this should be lambdaL -> inf, and its effect is to reduce z_real to 1/sigma_e and null z_imag,
        # we force this condition when lambdaL == 0
        lambdaL_zero_and_sigma_e_nonzero = not(lambdaL_nonzero) & sigma_e_nonzero;
        z_real(lambdaL_zero_and_sigma_e_nonzero) = (1.0 ./ sigma_e)(lambdaL_zero_and_sigma_e_nonzero);
        z_imag(lambdaL_zero_and_sigma_e_nonzero) = 0.0;
        # Special case: superconductor with only Cooper pairs
        # for sigma_e = 0 but lambdaL finite, z_real is null and z_imag reduces to 1j*omega*mu*lambdaL^2;
        # however this condition is already covered by in the calculation above, with no numerical issues
        % (den is 1, sigma_e = 0 means z_real is null)
    end
    
    % circulant tensor for the current frequency

    % no need to compute FFT every time, we already calculated it once for all.
    % FFT is a linear operator, and as the frequency-dependent
    % part is not in the circulant tensor, we can just multiply on by 'ko^2'
    #fN_all = fN_all2*(ko^2);
    %fN_all = fft_operator(fN_all);
    %st_sparse_precon = st_sparse_precon2 * (ko^2);
    % we can just multiply by 1j*omega*mu
    fN_all = (1j*omega*mu)*fN_all2;
    st_sparse_precon = (1j*omega*mu)*st_sparse_precon2;
    
    sim_CPU_lse(freq_no,1,1)=toc(tinisim); % CPU time for FFT + prep data
 
    % prepare the preconditioner
    tinisim = tic;
    lse_sparse_precon_prepare(dx,freq,z_real,z_imag,idxS,st_sparse_precon,nodeid_4_grnd,nodeid_4_injectcurr,Ae);
    sim_CPU_lse(freq_no,port_no,2)=toc(tinisim); % CPU time for sparse_precon
        
    for port_no=1:num_ports
        disp(['Solving for port # ',num2str(port_no), ' ...'])
        
        % ------------------------------------------------------------------------
        %     Solve Linear System of Equations Iteratively
        % -------------------------------------------------------------------------
        
        tinisim = tic;
        % Solve the system iteratively
        % Define the handle for matvect (remark: all other parameters beyond 'J'
        % are assigned at *this* time and are not modified any more later on when the function is called)
        fACPU   = @(J)lse_matvect_mult(J, fN_all, Ae, z_real, z_imag, dx, freq, idxS5, nodeid_4_grnd, nodeid_4_injectcurr);
        % Define the handle for the preconditioner multiplication (remark: all other parameters beyond 'JOut_full_in'
        % are assigned at *this* time and not modified any more later on when the function is called)
        fPCPU = @(JOut_full_in)lse_sparse_precon_multiply(JOut_full_in, Ae, nodeid_4_grnd, nodeid_4_injectcurr, prectol);
        tini = tic;
        disp(['Iterative solution started ... '])
        %[rhs_vect_sparse_precon]=lse_sparse_precon_multiply(rhs_vect(:,port_no),Ae,nodeid_4_grnd,nodeid_4_injectcurr);
        if (strcmp(fl_precon_type, 'schur_gmres') == 0)
            [x, flag, relres, iter, resvec] = pgmres(@(J)fACPU(J), rhs_vect(:,port_no), inner_it, tol, outer_it, @(JOut_full_in)fPCPU(JOut_full_in) );
        else
            [x, flag, relres, iter, resvec] = fpgmres(@(J)fACPU(J), rhs_vect(:,port_no), inner_it, tol, outer_it, @(JOut_full_in)fPCPU(JOut_full_in) );
        end

        tend = toc(tini);
        disp(['Total time for iterative solution ::: ' ,num2str(tend)]);
        disp(['Done... Iterative solution'])
        sim_CPU_lse(freq_no,port_no,3)=toc(tinisim); % CPU time for iterative solver
        if (abs(freq_curr_plot-freq)<1e-12 && port_no == 1)
            x_backup = x;
        end
        % ------------------------------------------------------------------------
        %     Compute the Currents on Port Nodes and the Column in Ymat
        % -------------------------------------------------------------------------

        currs_port_yparams=zeros(num_ports,1);
        for kk=1:num_ports
            currs_port_yparams(kk,1)=sum(rhs_vect(:,kk).*x);
        end
        
        Y_mat(:,port_no,freq_no)=currs_port_yparams(:,1);
        
        % remove the following later on! not needed!
        % compute column w/ alternative way - just for double checking
        %[currs_port_yparams2] = lse_compute_Y_mat_column_alternative(num_ports,Ae,Ae_only_leaving,Ae_only_entering_bndry,x,nodeid_lft,currs_port_yparams);
  
        
        disp(['Done... Solving for port # ',num2str(port_no)])
    end
    disp(['Done... Simulation for frequency : ',num2str(freq),' freq pnt: ',num2str(freq_no), ' / ', num2str(num_freq)])
    
    Z_mat(:,:,freq_no)=inv(squeeze(Y_mat(:,:,freq_no)));
    R_jL_mat(:,:,freq_no)=abs(real(squeeze(Z_mat(:,:,freq_no))))+sqrt(-1)*abs(imag(Z_mat(:,:,freq_no))/(2*pi*freq));
end

disp('-----------------------------------------------------')
disp(['Done... Solving LSEs'])
disp('-----------------------------------------------------')

for freq_no=1:num_freq
    disp(['R+jL matrix for frequency = ',num2str(freq_all(freq_no))])
    for kk=1:num_ports
        disp([num2str(R_jL_mat(kk,:,freq_no))])
    end
end

% ------------------------------------------------------------------------
%                         Storing Data
% -------------------------------------------------------------------------


disp('-----------------------------------------------------')
disp(['Saving Data...'])

% R+jL matrices
disp(['  Saving R+jL matrices'])
save(['Results', filesep, fileinname, '-', 'data_R_jL_mat.mat'], 'num_freq', 'num_ports','freq_all','R_jL_mat');

% CPU timings
disp(['  Saving CPU tinings'])
save(['Results', filesep, fileinname, '-', 'data_CPU_timings.mat'], 'num_freq', 'num_ports','sim_CPU_pre','sim_CPU_lse');

% Current plots
disp(['  Saving voxel currents for plotting'])
save(['Results', filesep, fileinname, '-', 'data_curr_plot.mat'], 'x', 'Ae_only_leaving','Ae_only_entering_bndry','Mc','dx','plot_option');

disp(['Done... Saving data'])
disp('-----------------------------------------------------')

if(plot_currents_post_proc == 1)
    %% ------------------------------------------------------------------------
    %                         Visualization
    % -------------------------------------------------------------------------

    disp('-----------------------------------------------------')
    disp(['Plotting Current Distribution...'])

    %close all

    % select plotting option - check the subroutines below for more options
    % option 1-> total currents on 3D structure, - no cut selection required
    % option 2-> current coefficients on the voxels - select plane and cut
    % option 3-> currents on the nodes via imagesc - select plane and cut
    % option 4-> currents on the nodes via quiver - select plane and cut
    % option 5-> currents on the structure w/directions via quiver3 - no cut selection required
    % voxels (on a selected cut), 3-> currents on nodes w/scalar values (on a selected cut) ,
    % plot_option=1;
    
    % generate domain 3D grid
    bbox_min=[0 0 0]; % minimum coordinates of bounding box
    bbox_max=[dx*L dx*M dx*N]; % max coordinates of bbox
    [r] = generategridfrombbox(dx, [bbox_min(1) bbox_max(1)],[bbox_min(2) bbox_max(2)],[bbox_min(3) bbox_max(3)],0);

    % set x_backup to x
    x = x_backup;

    % if any of plot option 2,3,4 is selected, define plane and cut
    slct_plane='xy'; %'xz'; 'yz';
    if (plot_option == 2 || plot_option == 3)
        % 1) use the following for plot option 2 and 3
        slct_cut=round(N/2);% round(M/2); round(L/2);
    elseif (plot_option == 4)
        % 2) use the following for plot option 4 - we need coordinate of the cut
        slct_cut=squeeze(r(1,1,N,3)); % z-coordinate of cut % squeeze(r(round(L/2),1,1,1)); % x-coordinate of cut; squeeze(r(1,round(M/2),1,2)); % y-coordinate of cut
    end

    if (plot_option == 2)
        % sort current coefficients on voxels
        [Jx_currs_grid,Jy_currs_grid,Jz_currs_grid,J2d_currs_grid,J3d_currs_grid,cmin,cmax]=post_obtain_curr_coefs_on_grid(x,Mc);
    elseif (plot_option > 2)
        % obtain currents on nodes
        [nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned]=post_obtain_currs_on_nodes(x,Ae_only_leaving,Ae_only_entering_bndry,r,Mc,dx);
    end

    switch plot_option
        case 1
            % 1) Plot currents on structure
            % Plot total currents as one scalar on each voxel
            plot_currs_on_3D_structure(x,Ae_only_leaving,r,Mc,dx)
        case 2
            % 2) Plot current coefficients obtained via iterative solution
            plot_curr_coefs_on_grid(slct_plane,slct_cut,r,Jx_currs_grid,Jy_currs_grid,Jz_currs_grid,J2d_currs_grid,J3d_currs_grid,cmin,cmax);
        case 3
            % 3) Plot currents with scalar values via imagesc
            plot_curr_on_nodes(slct_plane,slct_cut,dx,nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned)
        case 4
            % 4) Plot currents on cuts w/ directions via quiver
            plot_curr_on_nodes_quiver(slct_plane,slct_cut,nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned)
        case 5
            % 5) Plot currents on the structure w/directions via quiver3
            plot_curr_on_nodes_quiver3(nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned)
        otherwise
            disp('No current plotting!')
    end

    disp(['Done... Plotting Current Distribution'])
    disp('-----------------------------------------------------')

end
