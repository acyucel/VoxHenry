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
% straight conductor
%

%
% These are the user parameters. Modify this section to generate different models
%

% voxel size (deltax)
Res = 0.5e-6;
% inputs for generating conductors with specified lengths and widths of arms
num_conds = 2; % number of conductors
num_ports = 2; % number of ports
len_cond=30.0e-6; % length of conductors
width_cond=10.0e-6; % width of conductor
height_cond=10.0e-6; % height of conductor
dist_btw_conds=20.0e-6; % distance between centers of conductors

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
%freq = [1e10] ; % frequency
num_freq = length(freq);
freq_curr_plot=2.5e9; % frequency for plotting currents

er = 0;  % epsilon_r of conductors
%se=[5.8e7 3.77e7]; % conductivity of conductors
se=5.8e7; % conductivity of conductors
fl_check_domain=0; % set to 1 for only plotting the structure (no simulation)
fl_check_geo=0; % set to 1 for only plotting the domain (no simulation)
fl_check_ports=0; % set to 1 for only plotting the port nodes (no simulation)
plot_option=1; % see the options of plotting in Visualization part


% output file name
fileoutname = sprintf('Input_files%sstraight_cond%d_len%.1fu_wid%.1fu_dist%.1fu.vhr', filesep, num_conds, len_cond*1e6, width_cond*1e6, dist_btw_conds*1e6)

cen_cond1=[len_cond/2 width_cond/2 height_cond/2];
Cnt = [cen_cond1;]; % centers of conductors
for kk=1:num_conds-1
    for ll=1:1
        Cnt=[Cnt; [Cnt(ll,1) Cnt(ll,2)+kk*dist_btw_conds Cnt(ll,3)]];
    end
end

Dims_tmp = [len_cond width_cond height_cond;]; % dimensions of conductors(L(x),W(y),H(z))
Orients_tmp=['x';]; % orientations of conductors
Dims=[]; Orients=num2str([]);
for kk=1:num_conds
    Dims=[Dims;Dims_tmp];
    Orients=[Orients;Orients_tmp];
end

% -------------------------------------------------------------------------
%                         Initialize stuff
% -------------------------------------------------------------------------


disp(['Model type: straight conductors with width ', num2str(width_cond), ' m, height ', num2str(height_cond), ' m and length ', num2str(len_cond), ' m']);
disp(['            number of conductors: ', num2str(num_conds), ', distance between conductors ', num2str(dist_btw_conds)]);
disp(['            number of ports: ', num2str(num_ports)]);
disp('------------------------------------------------------------------')

firstline = sprintf('* Model type: straight conductors, width %g m, height %g m, length %g m\n*             distance %g, number of conductors %d, ports %d\n', width_cond, height_cond, len_cond, dist_btw_conds, num_conds, num_ports);

if (issorted(freq) == 0) % not sorted
    freq=sort(freq)
end
freq_all = freq;
freq = freq(1); % currently do everything for the lowest freq

EMconstants

bbox_min=[0 0 0]; % minimum coordinates of bounding box (bbox) - set to positive reals if possible
bbox_max=[len_cond (num_conds-1)*dist_btw_conds+width_cond height_cond]; % max coordinates of bbox

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
        pnt_rght{1}(dum,1:3)=[len_cond (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res)]; % points on which excitation defined
        pnt_lft{1}(dum,1:3)=[0 (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res)]; % points on which ground defined
        dum=dum+1;
    end
end

% defining the nodes in remaining ports
for kk=2:num_ports
    pnt_lft{kk}(:,2)=pnt_lft{1}(:,2) + dist_btw_conds*(kk-1);
    pnt_lft{kk}(:,[1 3])=pnt_lft{1}(:,[1 3]);
    pnt_rght{kk}(:,2)=pnt_rght{1}(:,2) + dist_btw_conds*(kk-1);
    pnt_rght{kk}(:,[1 3])=pnt_rght{1}(:,[1 3]);
end

% defining nodes connected ground if conductors without ports exist; if
% there is no, then leave as a empty array.

%pnt_well_cond=[pnt_lft{1}(1,1) pnt_lft{1}(1,2) + dist_btw_conds pnt_lft{1}(1,3);];
pnt_well_cond=[];

% -------------------------------------------------------------------------
%                   Define domain and constitutive parameters
% -------------------------------------------------------------------------

% generate domain 3D grid
[r] = generategridfrombbox(Res,[bbox_min(1) bbox_max(1)],[bbox_min(2) bbox_max(2)],[bbox_min(3) bbox_max(3)],fl_check_domain);

% assign constitutive parameters
[idx,epsilon_r,sigma_e,grid_intcon] = intcon_constparams(r,Res,Cnt,Dims,Orients,er,se,fl_check_geo);


% -------------------------------------------------------------------------
%           Output file
% -------------------------------------------------------------------------

pre_output_file

