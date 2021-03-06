clc; close all; clear all; format long e;

disp('------------------------------------------------------------------')
disp('                VoxHenry input file generator')
disp('------------------------------------------------------------------')

% -------------------------------------------------------------------------
%                  Add the Current Path to Workspace
% -------------------------------------------------------------------------

pre_define_the_path_for_folders

% -------------------------------------------------------------------------
%                  Inputs for the Structure
% -------------------------------------------------------------------------

%   
% square coil
%

%
% These are the user parameters. Modify this section to generate different models
%

% voxel size (deltax)
Res = 0.25e-6;

% inputs for generating conductors with specified lengths and widths of arms
num_conds = 1; % number of conductors
num_ports = 1; % number of ports

len_arm=100.0e-6; % the length of coil's arm 
width_cond=5.0e-6; % width of conductor
height_cond=5.0e-6; % height of conductor
spc_for_exc=2.0e-6; % spacing for excitation

% 
%  End of user parameters
% 

% -------------------------------------------------------------------------
%                  Inputs for Simulation
% -------------------------------------------------------------------------
freq = [1e0 2.5e0 5e0 7.5e0 1e1 2.5e1 5e1 7.5e1 1e2 2.5e2 5e2 7.5e2 ...
     1e3 2.5e3 5e3 7.5e3 1e4 2.5e4 5e4 7.5e4 1e5 2.5e5 5e5 7.5e5 ...
     1e6 2.5e6 5e6 7.5e6 1e7 2.5e7 5e7 7.5e7 1e8 2.5e8 5e8 7.5e8 ...
     1e9 2.5e9 5e9 7.5e9 1e10] ; % frequency

num_freq = length(freq);
freq_curr_plot=2.5e9; % frequency for plotting currents

er = 0;  % epsilon_r of conductors
se=5.8e7; % conductivity of conductors
lL=0.0; % London penetration depth (zero for non superconductive)
fl_check_domain=0; % set to 1 for only plotting the structure (no simulation)
fl_check_geo=0; % set to 1 for only plotting the domain (no simulation)
fl_check_ports=0; % set to 1 for only plotting the port nodes (no simulation)
plot_option=1; % see the options of plotting in Visualization part


% output file name
fileoutname = sprintf('Input_files%ssquare_coil_len%.1fu_wid%.1fu_heig%.1f_dist%.1fu.vhr', filesep, len_arm*1e6, width_cond*1e6, height_cond*1e6, spc_for_exc*1e6)

cen_cond1=[(len_arm/2-spc_for_exc/2)/2 height_cond/2 width_cond/2];
cen_cond2=[(len_arm+(len_arm/2+spc_for_exc/2))/2 height_cond/2 width_cond/2];
cen_cond3=[len_arm-(width_cond/2) height_cond/2 len_arm/2];
cen_cond4=[len_arm/2 height_cond/2 len_arm-(width_cond/2)];
cen_cond5=[width_cond/2 height_cond/2 len_arm/2];
Cnt = [cen_cond1; cen_cond2; cen_cond3; cen_cond4; cen_cond5;]; % centers of conductors
Dims_tmp1 = [len_arm height_cond width_cond;]; % dimensions of conductors(L(x),W(y),H(z))
Dims_tmp2 = [(len_arm-spc_for_exc)/2 height_cond width_cond;]; % dimensions of conductors(L(x),W(y),H(z))
Dims=[Dims_tmp2; Dims_tmp2; Dims_tmp1;Dims_tmp1;Dims_tmp1;];
Orients=['x';'x';'z';'x';'z';]; % orientations of conductors

% -------------------------------------------------------------------------
%                         Initialize stuff
% -------------------------------------------------------------------------


disp(['Model type: square coil with coil arm length ', num2str(len_arm), ' m, conductor width ', num2str(width_cond), ' m,']);
disp(['            conductor height ', num2str(height_cond), ' m, spacing ', num2str(spc_for_exc), ' m,']);
disp('------------------------------------------------------------------')

firstline = sprintf('* Model type: square coil with arm length %g m, width %g m, height %g m, spacing %g m\n', len_arm, width_cond, height_cond, spc_for_exc);

if (issorted(freq) == 0) % not sorted
    freq=sort(freq)
end
freq_all = freq;
freq = freq(1); % currently do everything for the lowest freq

EMconstants


% -------------------------------------------------------------------------
%                  Input for Computational Domain
% -------------------------------------------------------------------------
% At the end of this part, we only need bbox_min(3) and bbox_max(3) vectors
% define computational domain or bounding box enclosing the structure
bbox_min=[0 0 0]; % minimum coordinates of bounding box (bbox) - set to positive reals if possible
bbox_max=[len_arm width_cond len_arm]; % max coordinates of bbox

% generate domain 3D grid
[r] = generategridfrombbox(Res,[bbox_min(1) bbox_max(1)],[bbox_min(2) bbox_max(2)],[bbox_min(3) bbox_max(3)],fl_check_domain);

% assign constitutive parameters
[idx,epsilon_r,sigma_e,lambdaL,grid_intcon] = intcon_constparams(r,Res,Cnt,Dims,Orients,er,se,lL,fl_check_geo);

% -------------------------------------------------------------------------
%                  Input for Ports
% -------------------------------------------------------------------------
% At the end of this part, we need structures pnt_lft{xx} and pnt_rght{xx} which
% contains the coordinates of nodes on both sides of xxth port

% defining the nodes in first port
pnt_lft=cell(num_ports,1);
pnt_rght=cell(num_ports,1);
pnt_lft{1}=zeros(round(width_cond/Res)*round(height_cond/Res),3);
pnt_rght{1}=zeros(round(width_cond/Res)*round(height_cond/Res),3);
dum=1;
for kk=1:round(width_cond/Res)
    for ll=1:round(height_cond/Res)
        pnt_rght{1}(dum,1:3)=[len_arm/2-spc_for_exc/2 (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res)]; % points on which excitation defined
        pnt_lft{1}(dum,1:3)=[len_arm/2+spc_for_exc/2 (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res)]; % points on which ground defined
        dum=dum+1;
    end
end

% defining nodes connected ground if conductors without ports exist; if
% there is no, then leave as a empty array.

%pnt_well_cond=[pnt_lft{1}(1,1) pnt_lft{1}(1,2) + dist_btw_conds pnt_lft{1}(1,3);];
pnt_well_cond=[];

% -------------------------------------------------------------------------
%           Output file
% -------------------------------------------------------------------------

pre_output_file
