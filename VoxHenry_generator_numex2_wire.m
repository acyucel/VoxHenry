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
% wire
%

%
% These are the user parameters. Modify this section to generate different models
%

% voxel size (deltax)
Res = 1.0e-6; 
% inputs for generating conductors with specified lengths and radii 
num_ports = 1; % number of ports
len_cond=50.0e-6; % length of conductors
dia_cond=10.0e-6; % diameter of conductor
%len_cond=2.0e-6; % length of conductors
%dia_cond=4.0e-6; % diameter of conductor

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
se=5.8e7; % conductivity of conductors
lL=0.0; % London penetration depth (zero for non superconductive)
fl_check_domain=0; % set to 1 for only plotting the structure (no simulation)
fl_check_geo=0; % set to 1 for only plotting the domain (no simulation)
fl_check_ports=0; % set to 1 for only plotting the port nodes (no simulation)
plot_option=1; % see the options of plotting in Visualization part


% output file name
fileoutname = sprintf('Input_files%swire_len%.1fu_dia%.1fu.vhr', filesep, len_cond*1e6, dia_cond*1e6)

rad_cond=dia_cond/2; % radius of conductor
cen_cond1=[len_cond/2 dia_cond/2 dia_cond/2];
Cnt = [cen_cond1;]; % centers of conductors

Dims_tmp = [len_cond dia_cond dia_cond;]; % dimensions of conductors(L(x),W(y),H(z))
Orients_tmp=['x';]; % orientations of conductors
Dims=[]; Orients=[];
Dims=[Dims;Dims_tmp];
Orients=[num2str(Orients);Orients_tmp];


% -------------------------------------------------------------------------
%                         Initialize stuff
% -------------------------------------------------------------------------


disp(['Model type: cylindrical wire with diameter ', num2str(dia_cond), ' m and length ', num2str(len_cond), ' m']);
disp('------------------------------------------------------------------')

firstline = sprintf('* Model type: cylindrical wire with radius %g m and length %g m\n', dia_cond, len_cond);

if (issorted(freq) == 0) % not sorted
    freq=sort(freq)
end
freq_all = freq;
freq = freq(1); % currently do everything for the lowest freq

EMconstants


% -------------------------------------------------------------------------
%                   Define domain and constitutive parameters
% -------------------------------------------------------------------------

% define bounding box enclosing the structure
bbox_min=[0 0 0]; % minimum coordinates of bounding box (bbox) - set to positive reals if possible
bbox_max=[len_cond dia_cond dia_cond]; % max coordinates of bbox

% generate domain 3D grid
[r] = generategridfrombbox(Res,[bbox_min(1) bbox_max(1)],[bbox_min(2) bbox_max(2)],[bbox_min(3) bbox_max(3)],fl_check_domain);

% assign constitutive parameters
[idx,epsilon_r,sigma_e,lambda_L,grid_intcon] = intcon_constparams(r,Res,Cnt,Dims,Orients,er,se,lL,fl_check_geo);

% removing the elements outside of circle of each cross-sectional area
anchor_pnt_y_z=[rad_cond rad_cond];
for kk=1:size(grid_intcon,1)
    anchor_pnt=[squeeze(r(kk,1,1,1)) anchor_pnt_y_z]';
    for ll=1:size(grid_intcon,2)
        for mm=1:size(grid_intcon,3)
            if (norm(squeeze(grid_intcon(kk,ll,mm,1:3))-anchor_pnt) > rad_cond)
                epsilon_r(kk,ll,mm)=1;
                sigma_e(kk,ll,mm)=0;
                if any(lL)
                    lambda_L(kk,ll,mm)=0;
                end
                grid_intcon(kk,ll,mm,1:3)=0;
            end
                
        end
    end
end
        
if (fl_check_domain == 1)
    plot_boxes_of_grid(grid_intcon,Res);
end

if (fl_check_domain == 1 || fl_check_geo == 1); return; end;


% -------------------------------------------------------------------------
%                  Input for Ports
% -------------------------------------------------------------------------
% At the end of this part, we need structures pnt_lft{xx} and pnt_rght{xx} which
% contains the coordinates of nodes on both sides of xxth port

% defining the nodes in first port
pnt_lft=cell(num_ports,1);
pnt_rght=cell(num_ports,1);
dum=1;
for kk=1:round(dia_cond/Res)
    for ll=1:round(dia_cond/Res)
        if (abs(grid_intcon(1,kk,ll,1)) > 1e-12 || ...
                abs(grid_intcon(1,kk,ll,2)) > 1e-12 || ...
                abs(grid_intcon(1,kk,ll,3)) > 1e-12)
            tmp_coor=squeeze(grid_intcon(1,kk,ll,1:3));
            pnt_rght{1}(dum,1:3)=[tmp_coor(1)-Res/2 tmp_coor(2) tmp_coor(3)]; % points on which excitation defined
            pnt_lft{1}(dum,1:3)=[tmp_coor(1)-Res/2+len_cond tmp_coor(2) tmp_coor(3)]; % points on which ground defined
            dum=dum+1;
        end
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
