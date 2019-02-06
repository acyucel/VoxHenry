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
% circular coil
%

%
% These are the user parameters. Modify this section to generate different models
%

% voxel size (deltax)
Res = 1.0e-6; 
% inputs for generating circular coil with specified radii
num_ports = 1; % number of ports

rad_loop = 150.0e-6; % radius of coil
rad_wire = 5.0e-6; % % radius of tube

% 
%  End of user parameters
% 

% -------------------------------------------------------------------------
%                  Inputs for Simulation
% -------------------------------------------------------------------------
freq = [1e0 2.5e0 5e0 7.5e0 1e1 2.5e1 5e1 7.5e1 1e2 2.5e2 5e2 7.5e2 ...
    1e3 2.5e3 5e3 7.5e3 1e4 2.5e4 5e4 7.5e4 1e5 2.5e5 5e5 7.5e5 ...
    1e6 2.5e6 5e6 7.5e6 1e7 2.5e7 5e7 7.5e7 1e8 2.5e8 5e8 7.5e8 ...
    1e9 2.5e9 5e9 7.5e9 1e10 2.5e10 5e10 7.5e10 1e11 2.5e11 5e11 7.5e11 1e12] ; % frequency

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
fileoutname = sprintf('Input_files%scircular_coil_loop%.1fu_rad%.1fu.vhr', filesep, rad_loop*1e6, rad_wire*1e6)

cen_cond1=[(rad_wire+rad_loop) (rad_wire+rad_loop) rad_wire];
Cnt = [cen_cond1;]; % centers of conductors
Dims_tmp1 = [2*(rad_wire+rad_loop) 2*(rad_wire+rad_loop) 2*rad_wire]; % dimensions of conductors(L(x),W(y),H(z))
Dims=[Dims_tmp1;];
Orients=['x';]; % orientations of conductors

% -------------------------------------------------------------------------
%                         Initialize stuff
% -------------------------------------------------------------------------


disp(['Model type: circular coil with coil radius ', num2str(rad_loop), ' m and wire radius ', num2str(rad_wire), ' m']);
disp('------------------------------------------------------------------')

firstline = sprintf('* Model type: circular coil with coil radius %g m and wire radius %g m\n', rad_loop, rad_wire);

if (issorted(freq) == 0) % not sorted
    freq=sort(freq)
end
freq_all = freq;
freq = freq(1); % currently do everything for the lowest freq

EMconstants

% -------------------------------------------------------------------------
%                   Define domain and constitutive parameters
% -------------------------------------------------------------------------

% define computational domain bounding box enclosing the structure
bbox_min=[0 0 0]; % minimum coordinates of bounding box (bbox) - set to positive reals if possible
bbox_max=[2*(rad_wire+rad_loop)+1e-12 2*(rad_wire+rad_loop)+1e-12 2*rad_wire+1e-12]; % max coordinates of bbox

% generate domain 3D grid
[r] = generategridfrombbox(Res,[bbox_min(1) bbox_max(1)],[bbox_min(2) bbox_max(2)],[bbox_min(3) bbox_max(3)],fl_check_domain);

% assign constitutive parameters
[idx,epsilon_r,sigma_e,lambdaL,grid_intcon] = intcon_constparams(r,Res,Cnt,Dims,Orients,er,se,lL,fl_check_geo);

% removing the elements outside of circular loop (torus - equation)
shft_org_torrus=[(rad_wire+rad_loop) (rad_wire+rad_loop) rad_wire];
for kk=1:size(grid_intcon,1)
    for ll=1:size(grid_intcon,2)
        for mm=1:size(grid_intcon,3)
            tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3))'-shft_org_torrus;
            equa_tmp=(rad_loop-sqrt(tmp_coor(1)^2+tmp_coor(2)^2))^2+(tmp_coor(3)^2);
            if (equa_tmp >= rad_wire^2)
                epsilon_r(kk,ll,mm)=1;
                sigma_e(kk,ll,mm)=0;
                if any(lL)
                    lambdaL(kk,ll,mm)=0;
                end
                grid_intcon(kk,ll,mm,1:3)=0;
            end
            tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3))';
            if (abs(tmp_coor(1)-((rad_wire+rad_loop)+0.5*Res)) < 1e-12 && ...
                    tmp_coor(2) < ((rad_wire+rad_loop)+0.5*Res))
                epsilon_r(kk,ll,mm)=1;
                sigma_e(kk,ll,mm)=0;
                if any(lL)
                    lambdaL(kk,ll,mm)=0;
                end
                grid_intcon(kk,ll,mm,1:3)=0;
            end
        end
    end
end

if (fl_check_geo == 1)
    tic
    plot_boxes_of_grid(grid_intcon,Res);
    disp(['Time for visualizing the structure ::: ',num2str(toc)])
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
ind_box_exc=round((((rad_wire+rad_loop)-Res/2)-Res/2)/(Res))+1;
dum=1;
if isempty(lambdaL)
    for kk=ind_box_exc:ind_box_exc
        for ll=1:size(grid_intcon,2)
            for mm=1:size(grid_intcon,3)
                if (sigma_e(kk,ll,mm) > 0 && squeeze(grid_intcon(kk,ll,mm,2)) < (rad_wire+rad_loop)) % conductor
                    tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3));
                    pnt_lft{1}(dum,1:3)=[tmp_coor(1)+Res/2 tmp_coor(2) tmp_coor(3)];
                    dum=dum+1;
                end
            end
        end
    end
else
    for kk=ind_box_exc:ind_box_exc
        for ll=1:size(grid_intcon,2)
            for mm=1:size(grid_intcon,3)
                if ( (sigma_e(kk,ll,mm) > 0 || lambdaL(kk,ll,mm) > 0) && squeeze(grid_intcon(kk,ll,mm,2)) < (rad_wire+rad_loop)) % conductor
                    tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3));
                    pnt_lft{1}(dum,1:3)=[tmp_coor(1)+Res/2 tmp_coor(2) tmp_coor(3)];
                    dum=dum+1;
                end
            end
        end
    end
end

ind_box_grnd=round((((rad_wire+rad_loop)+3*Res/2)-Res/2)/(Res))+1;
dum=1;
if isempty(lambdaL)
    for kk=ind_box_grnd:ind_box_grnd
        for ll=1:size(grid_intcon,2)
            for mm=1:size(grid_intcon,3)
                if (sigma_e(kk,ll,mm) > 0 && squeeze(grid_intcon(kk,ll,mm,2)) < (rad_wire+rad_loop) ) % conductor
                    tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3));
                    pnt_rght{1}(dum,1:3)=[tmp_coor(1)-Res/2 tmp_coor(2) tmp_coor(3)];
                    dum=dum+1;
                end
            end
        end
    end
else
    for kk=ind_box_grnd:ind_box_grnd
        for ll=1:size(grid_intcon,2)
            for mm=1:size(grid_intcon,3)
                if ((sigma_e(kk,ll,mm) > 0 || lambdaL(kk,ll,mm) > 0) && squeeze(grid_intcon(kk,ll,mm,2)) < (rad_wire+rad_loop) ) % conductor
                    tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3));
                    pnt_rght{1}(dum,1:3)=[tmp_coor(1)-Res/2 tmp_coor(2) tmp_coor(3)];
                    dum=dum+1;
                end
            end
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
