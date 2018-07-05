%%%%%%%%%%%%%%% Several Structures for EM Analysis %%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Square Coil %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Square coil with specified arm width and length on xz plane (extruded along y direction)

freq = 3e9 %3e12; % frequency
% inputs for generating square coil with specified lengths and widths of arms
len_arm=50.0e-6;
width_arm=5.0e-6;
height_arm=5.0e-6; 

% auto square coil generation
% Resolution
Res = 5.0e-6;
bbox_size=[len_arm height_arm len_arm]; % domain size
cen_cond1=[(len_arm/2-Res)/2 height_arm/2 width_arm/2];
cen_cond2=[(len_arm+(len_arm/2+Res))/2 height_arm/2 width_arm/2];
cen_cond3=[len_arm-(width_arm/2) height_arm/2 len_arm/2];
cen_cond4=[len_arm/2 height_arm/2 len_arm-(width_arm/2)];
cen_cond5=[width_arm/2 height_arm/2 len_arm/2];
Cnt = [cen_cond1; cen_cond2; cen_cond3; cen_cond4; cen_cond5;]; % centers of conductors
Dims_tmp1 = [len_arm height_arm width_arm;]; % dimensions of conductors(L(x),W(y),H(z))
Dims_tmp2 = [(len_arm-Res)/2 height_arm width_arm;]; % dimensions of conductors(L(x),W(y),H(z))
Dims=[Dims_tmp2; Dims_tmp2; Dims_tmp1;Dims_tmp1;Dims_tmp1;];
Orients=['x';'x';'z';'x';'z';]; % orientations of conductors
er = 0;  % epsilon_r of interconnect
se=5.8e7; % conductivity of interconnect


% define excitation and ground nodes in the spacing between arms
dum=1;
pnt_exc=zeros(round(width_arm/Res)*round(height_arm/Res),3);
pnt_grnd=zeros(round(width_arm/Res)*round(height_arm/Res),3);
for kk=1:round(width_arm/Res)
    for ll=1:round(height_arm/Res)
        pnt_exc(dum,1:3)=[len_arm/2+Res (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res)]; % points on which excitation defined
        pnt_grnd(dum,1:3)=[len_arm/2-Res (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res)]; % points on which ground defined
        dum=dum+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Circular Wire %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2) Circular wire with specified radius and length

len_wire=50.0e-6;
rad_wire=10.0e-6;

% auto wire generation
% Resolution
Res =1.0e-6; %0.25e-6;% 0.5e-6; %
bbox_size=[len_wire 2*rad_wire 2*rad_wire]; % domain size
cen_cond1=[len_wire/2 rad_wire rad_wire];
Cnt = [cen_cond1;]; % centers of conductors
Dims_tmp1 = [len_wire 2*rad_wire 2*rad_wire;]; % dimensions of conductors(L(x),W(y),H(z))
Dims=[Dims_tmp1;];
Orients=['x';]; % orientations of conductors
er = 0;  % epsilon_r of interconnect
se=5.8e7; % conductivity of interconnect

% -------------------------------------------------------------------------
%                         Initialize stuff
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;
disp('-----------------------------------------------------')
disp('******** TERACIS: TERAhertz CIrcuit Simulator *******')
disp('-----------------------------------------------------')
disp(['Conductivity / rel. permittivity :: ',num2str(se),' / ',num2str(er)])
tmp_er=er - 1j*se/(eo*omega);
disp(['Relative epsr (Re/Im) :: ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
% Skin depth:
skin_depth=1/sqrt(pi*freq*4*pi*1e-7*se);
disp(['Discretization should be at order of (skin depth) ', num2str(skin_depth)])
% -------------------------------------------------------------------------
%                   Define domain and constitutive parameters
% -------------------------------------------------------------------------

% generate domain 3D grid
[r] = generategridfrombbox(Res,[0 bbox_size(1)+1e-12],[0 bbox_size(2)+1e-12],[0 bbox_size(3)+1e-12],1);

% assign constitutive parameters
[idx,epsilon_r,sigma_e,grid_intcon] = intcon_constparams(r,Res,Cnt,Dims,Orients,er,se,1);

% removing the elements outside of circle of each cross-sectional area
anchor_pnt_y_z=[rad_wire rad_wire];
for kk=1:size(grid_intcon,1)
    anchor_pnt=[squeeze(r(kk,1,1,1)) anchor_pnt_y_z]';
    for ll=1:size(grid_intcon,2)
        for mm=1:size(grid_intcon,3)
            if (norm(squeeze(grid_intcon(kk,ll,mm,1:3))-anchor_pnt) > rad_wire)
                epsilon_r(kk,ll,mm)=1;
                sigma_e(kk,ll,mm)=0;
                grid_intcon(kk,ll,mm,1:3)=0;
            end
                
        end
    end
end
        
tic
plot_boxes_of_grid(grid_intcon,Res);
disp(['Time for visualizing the structure ::: ',num2str(toc)])

% define excitation and ground nodes

dum=1;
for kk=1:1
    for ll=1:size(grid_intcon,2)
        for mm=1:size(grid_intcon,3)
            if (abs(grid_intcon(kk,ll,mm,1)) > 1e-12 || abs(grid_intcon(kk,ll,mm,2)) > 1e-12 || abs(grid_intcon(kk,ll,mm,3)) > 1e-12)
                tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3));
                pnt_exc(dum,1:3)=[tmp_coor(1)-Res/2 tmp_coor(2) tmp_coor(3)];
                dum=dum+1;
            end
        end
    end
end

dum=1;
for kk=size(grid_intcon,1):size(grid_intcon,1)
    for ll=1:size(grid_intcon,2)
        for mm=1:size(grid_intcon,3)
            if (abs(grid_intcon(kk,ll,mm,1)) > 1e-12 || abs(grid_intcon(kk,ll,mm,2)) > 1e-12 || abs(grid_intcon(kk,ll,mm,3)) > 1e-12)
                tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3));
                pnt_grnd(dum,1:3)=[tmp_coor(1)+Res/2 tmp_coor(2) tmp_coor(3)];
                dum=dum+1;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Circular Coil %%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = 3e9 %3e12; % frequency
% inputs for generating a cirular cross-section wire along x direction
rad_loop=50.0e-6;%50.0e-6;
rad_wire=5.0e-6;
fl_geo_check=0;
fl_method_Ae_gen=2;
fl_volt_source=1; % 1 for voltage source, 0 for current source
slct_decomp_sch='lu_decomp'; %'no_decomp','lu_decomp','ldlt_decomp','chol_decomp'

% auto circular loop generation
% Resolution
Res =1.0e-6; %0.5e-6; %0.25e-6;% 0.5e-6; %
bbox_size=[2*(rad_wire+rad_loop) 2*(rad_wire+rad_loop) 2*rad_wire]; % domain size
cen_cond1=[(rad_wire+rad_loop) (rad_wire+rad_loop) rad_wire];
Cnt = [cen_cond1;]; % centers of conductors
Dims_tmp1 = [bbox_size(1) bbox_size(2) bbox_size(3);]; % dimensions of conductors(L(x),W(y),H(z))
Dims=[Dims_tmp1;];
Orients=['x';]; % orientations of conductors
er = 0;  % epsilon_r of interconnect
se=5.8e7; % conductivity of interconnect

% -------------------------------------------------------------------------
%                         Initialize stuff
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;
disp('-----------------------------------------------------')
disp('******** TERACIS: TERAhertz CIrcuit Simulator *******')
disp('-----------------------------------------------------')
disp(['Conductivity / rel. permittivity :: ',num2str(se),' / ',num2str(er)])
tmp_er=er - 1j*se/(eo*omega);
disp(['Relative epsr (Re/Im) :: ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
% Skin depth:
skin_depth=1/sqrt(pi*freq*4*pi*1e-7*se);
disp(['Discretization should be at order of (skin depth) ', num2str(skin_depth)])
% -------------------------------------------------------------------------
%                   Define domain and constitutive parameters
% -------------------------------------------------------------------------

% generate domain 3D grid
[r] = generategridfrombbox(Res,[0 bbox_size(1)+1e-12],[0 bbox_size(2)+1e-12],[0 bbox_size(3)+1e-12],0);

% assign constitutive parameters
[idx,epsilon_r,sigma_e,grid_intcon] = intcon_constparams(r,Res,Cnt,Dims,Orients,er,se,0);

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
                grid_intcon(kk,ll,mm,1:3)=0;
            end
            tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3))';
            if (abs(tmp_coor(1)-((rad_wire+rad_loop)+0.5*Res)) < 1e-12 && ...
                    tmp_coor(2) < ((rad_wire+rad_loop)+0.5*Res))
                epsilon_r(kk,ll,mm)=1;
                sigma_e(kk,ll,mm)=0;
                grid_intcon(kk,ll,mm,1:3)=0;
            end
        end
    end
end
        
%tic
%plot_boxes_of_grid(grid_intcon,Res);
%disp(['Time for visualizing the structure ::: ',num2str(toc)])

% define excitation and ground nodes

ind_box_exc=round((((rad_wire+rad_loop)-Res/2)-Res/2)/(Res))+1;
dum=1;
for kk=ind_box_exc:ind_box_exc
    for ll=1:size(grid_intcon,2)
        for mm=1:size(grid_intcon,3)
            if (sigma_e(kk,ll,mm) > 0 && squeeze(grid_intcon(kk,ll,mm,2)) < (rad_wire+rad_loop)) % conductor
                tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3));
                pnt_exc(dum,1:3)=[tmp_coor(1)+Res/2 tmp_coor(2) tmp_coor(3)];
                dum=dum+1;
            end
        end
    end
end

ind_box_grnd=round((((rad_wire+rad_loop)+3*Res/2)-Res/2)/(Res))+1;
dum=1;
for kk=ind_box_grnd:ind_box_grnd
    for ll=1:size(grid_intcon,2)
        for mm=1:size(grid_intcon,3)
            if (sigma_e(kk,ll,mm) > 0 && squeeze(grid_intcon(kk,ll,mm,2)) < (rad_wire+rad_loop) ) % conductor
                tmp_coor=squeeze(grid_intcon(kk,ll,mm,1:3));
                pnt_grnd(dum,1:3)=[tmp_coor(1)-Res/2 tmp_coor(2) tmp_coor(3)];
                dum=dum+1;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%% Square Coil w/ ground plane %%%%%%%%%%%%%%%%%%%%%%%%

%freq = 3e9 %3e12; % frequency
% inputs for generating square coil with specified lengths and widths of arms
len_arm=20.0e-6;
width_arm=2.0e-6; %5.0e-6;
height_arm=2.0e-6;
% inputs for generating ground plane with specified length, width, and
% height
len_gp=40.0e-6;
width_gp=2.0e-6;
height_gp=2.0e-6;

dist_btw_coil_gp=3.0e-6;

fl_geo_check=0;
fl_method_Ae_gen=2;
fl_volt_source=1; % 1 for voltage source, 0 for current source
slct_decomp_sch='lu_decomp'; %'no_decomp','lu_decomp','ldlt_decomp','chol_decomp'

% auto square coil generation
% Resolution
Res = 0.25e-6; %0.5e-6; % resolution should be smaller than half_spc_feed
half_spc_feed=0.5e-6;
bbox_size=[len_gp height_arm+height_gp+dist_btw_coil_gp len_gp]; % domain size - w/assumption gp larger than or equal to coil
shft_vect=[len_gp/2-len_arm/2 0 len_gp/2-len_arm/2];

cen_cond1=[(len_arm/2-half_spc_feed)/2 height_arm/2 width_arm/2]+shft_vect;
cen_cond2=[(len_arm+(len_arm/2+half_spc_feed))/2 height_arm/2 width_arm/2]+shft_vect;
cen_cond3=[len_arm-(width_arm/2) height_arm/2 len_arm/2]+shft_vect;
cen_cond4=[len_arm/2 height_arm/2 len_arm-(width_arm/2)]+shft_vect;
cen_cond5=[width_arm/2 height_arm/2 len_arm/2]+shft_vect;
cen_cond6=[len_gp/2 height_arm+dist_btw_coil_gp+0.5*height_gp len_gp/2]; %gp
Cnt = [cen_cond1; cen_cond2; cen_cond3; cen_cond4; cen_cond5; cen_cond6;]; % centers of conductors

Dims_tmp1 = [len_arm height_arm width_arm;]; % dimensions of conductors(L(x),W(y),H(z))
Dims_tmp2 = [(len_arm-2*half_spc_feed)/2 height_arm width_arm;]; % dimensions of conductors(L(x),W(y),H(z))
Dims_tmp3 = [len_gp height_gp len_gp;]; % dimensions of conductors(L(x),W(y),H(z))
Dims=[Dims_tmp2; Dims_tmp2; Dims_tmp1;Dims_tmp1;Dims_tmp1;Dims_tmp3];
Orients=['x';'x';'z';'x';'z';'x']; % orientations of conductors

%Orients=['x';'x';'z';'x';'z';]; % orientations of conductors
er = 0;  % epsilon_r of interconnect
se=5.8e7; % conductivity of interconnect

% -------------------------------------------------------------------------
%                         Initialize stuff
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;
disp('-----------------------------------------------------')
disp('******** TERACIS: TERAhertz CIrcuit Simulator *******')
disp('-----------------------------------------------------')
disp(['Conductivity / rel. permittivity :: ',num2str(se),' / ',num2str(er)])
tmp_er=er - 1j*se/(eo*omega);
disp(['Relative epsr (Re/Im) :: ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
% Skin depth:
skin_depth=1/sqrt(pi*freq*4*pi*1e-7*se);
disp(['Discretization should be at order of (skin depth) ', num2str(skin_depth)])
% -------------------------------------------------------------------------
%                   Define domain and constitutive parameters
% -------------------------------------------------------------------------

% generate domain 3D grid
[r] = generategridfrombbox(Res,[0 bbox_size(1)+1e-12],[0 bbox_size(2)+1e-12],[0 bbox_size(3)+1e-12],0);

% assign constitutive parameters
[idx,epsilon_r,sigma_e,grid_intcon] = intcon_constparams(r,Res,Cnt,Dims,Orients,er,se,0);
%return
% define excitation and ground nodes in the spacing between arms for first

dum=1;
pnt_exc=zeros(round(width_arm/Res)*round(height_arm/Res),3);
pnt_grnd=zeros(round(width_arm/Res)*round(height_arm/Res),3);
for kk=1:round(width_arm/Res)
    for ll=1:round(height_arm/Res)
        pnt_exc(dum,1:3)=[len_arm/2+half_spc_feed+shft_vect(1) (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res)+shft_vect(3)]; % points on which excitation defined
        pnt_grnd(dum,1:3)=[len_arm/2-half_spc_feed+shft_vect(1) (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res)+shft_vect(3)]; % points on which ground defined
        dum=dum+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%% RF Coil Array %%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
%                  Define the Model to Simulate
% -------------------------------------------------------------------------

freq = 1e11 %3e12; % frequency
% inputs for generating N turn RF inductor with specified width and spacing
% of interconnects + max edge length of inductor
%slct_exc_coil=1; % select coil for excitation
n_turns=3;
n_arr_elems=2; % number of array elements
spc_elem=50e-6; % spacing between elements
len_arm=40.0e-6;
width_arm=3.0e-6;
height_arm=width_arm;
spc_arm=1.0e-6; % spacing btw conductors
len_feed=4.0e-6;
spc_dg=1.0e-6; % spacing for delta-gap ! attention: for coarse disc, check the spacing

fl_geo_check=0;
fl_method_Ae_gen=2;
fl_volt_source=0; % 1 for voltage source, 0 for current source
slct_decomp_sch='lu_decomp'; %'no_decomp','lu_decomp','ldlt_decomp','chol_decomp'

%iterative solver inputs
inner_it = 100; outer_it = 10; tol=1e-8; %tol=1e-8;

% auto square coil generation
% Resolution
Res = 0.5e-6; %0.5e-6; % resolution should be smaller than half_spc_feed
bbox_size_upper_bnd=[(n_arr_elems-1)*spc_elem+len_arm (n_arr_elems-1)*spc_elem+len_arm height_arm]; % domain upper limits
bbox_size_lower_bnd=[0 0 -len_feed]; % domain lower limits

grd_spc=spc_arm+width_arm;

cen_cond1=zeros(n_turns,3); % left arms
for kk=1:n_turns
    cen_cond1(kk,1:3)=[(grd_spc*(kk-1))+width_arm/2 len_arm/2 height_arm/2];
end

cen_cond2=zeros(n_turns,3); % upper arms
for kk=1:n_turns
    cen_cond2(kk,1:3)=[len_arm/2 (len_arm-width_arm/2)-(grd_spc*(kk-1)) height_arm/2];
end

cen_cond3=zeros(n_turns,3); % right arms
for kk=1:n_turns
    cen_cond3(kk,1:3)=[(len_arm-width_arm/2)-(grd_spc*(kk-1)) len_arm/2+grd_spc/2 height_arm/2];
end

cen_cond4=zeros(n_turns,3); % lower arms
for kk=1:n_turns
    cen_cond4(kk,1:3)=[len_arm/2+grd_spc/2 (grd_spc+width_arm/2)+(grd_spc*(kk-1)) height_arm/2];
end

cen_cond5=[width_arm/2 width_arm/2 -len_feed/2]; % feed z aligned 1
cen_cond6=[(n_turns*grd_spc)+width_arm/2 (n_turns*grd_spc)+width_arm/2 -len_feed/2]; % feed z aligned 1
cen_cond7=[(((n_turns*grd_spc)+width_arm)/2-spc_dg/2)/2 width_arm/2 -len_feed+width_arm/2]; % feed x aligned 1
cen_cond8=[((((n_turns*grd_spc)+width_arm)/2+spc_dg/2)+(((n_turns*grd_spc)+width_arm)))/2 width_arm/2 -len_feed+width_arm/2]; % feed x aligned 2
cen_cond9=[(n_turns*grd_spc)+width_arm/2 ((n_turns*grd_spc)+width_arm)/2 -len_feed+width_arm/2]; % feed y aligned 1

cen_tmp=[cen_cond1;cen_cond2;cen_cond3;cen_cond4;cen_cond5;cen_cond6;cen_cond7;cen_cond8;cen_cond9;];

shft_vects=zeros(n_arr_elems^2,3);
dum=1;
for kk=1:n_arr_elems % for elems along x
    for ll=1:n_arr_elems % for elems along y
        shft_vects(dum,1:3)=[(kk-1)*spc_elem (ll-1)*spc_elem 0];
        dum=dum+1;
    end
end

Cnt=[];
for kk=1:n_arr_elems^2
    Cnt = [Cnt;[cen_tmp(:,1)+shft_vects(kk,1) cen_tmp(:,2)+shft_vects(kk,2) cen_tmp(:,3)+shft_vects(kk,3)]]; % centers of conductors
end

Dims_tmp1=zeros(n_turns,3); % left & upper arms
Dims_tmp2=zeros(n_turns,3); % right arms

Orients_tmp1=[];Orients_tmp2=[];
for kk=1:n_turns
    Dims_tmp1(kk,1:3) = [len_arm-((kk-1)*2*grd_spc) width_arm height_arm;];
    Dims_tmp2(kk,1:3) = [len_arm-((kk-1)*2*grd_spc)-grd_spc width_arm height_arm;];
    if (len_arm-((kk-1)*2*grd_spc) < 0)
        disp('reduce the number of turns or increase the length of inductor')
        return
    end
    Orients_tmp1=[Orients_tmp1;'y';];
    Orients_tmp2=[Orients_tmp2;'x';];
end

Dims_tmp3=[len_feed width_arm height_arm;]; % feed z-aligned comp
Dims_tmp4=[((n_turns*grd_spc)+width_arm-spc_dg)/2 width_arm height_arm;]; % feed x-aligned comp
Dims_tmp5=[(n_turns*grd_spc)+width_arm width_arm height_arm;]; % feed y-aligned comp

% left, upper, right, lower, feed z1, feed z2;
Dims_all_tmp=[Dims_tmp1;Dims_tmp1;Dims_tmp2;Dims_tmp2;Dims_tmp3;Dims_tmp3;Dims_tmp4;Dims_tmp4;Dims_tmp5;];
Orients_all_tmp=[Orients_tmp1;Orients_tmp2;Orients_tmp1;Orients_tmp2;'z';'z';'x';'x';'y';]; % orientations of conductors

Dims=[];
Orients=[];
for kk=1:n_arr_elems^2
    Dims = [Dims;Dims_all_tmp];
    Orients = [Orients; Orients_all_tmp];
end

er = 0;  % epsilon_r of interconnect
se=5.8e7; % conductivity of interconnect

% -------------------------------------------------------------------------
%                         Initialize Stuff
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;
disp('-----------------------------------------------------')
disp('******** TERACIS: TERAhertz CIrcuit Simulator *******')
disp('-----------------------------------------------------')
disp(['Conductivity / rel. permittivity :: ',num2str(se),' / ',num2str(er)])
tmp_er=er - 1j*se/(eo*omega);
disp(['Relative epsr (Re/Im) :: ',num2str(real(tmp_er)),' / ',num2str(-imag(tmp_er))])
% Skin depth:
skin_depth=1/sqrt(pi*freq*4*pi*1e-7*se);
disp(['Discretization should be at order of (skin depth) ', num2str(skin_depth)])
% -------------------------------------------------------------------------
%                   Define Domain and Constitutive Parameters
% -------------------------------------------------------------------------

% generate domain 3D grid
[r] = generategridfrombbox(Res,[bbox_size_lower_bnd(1) bbox_size_upper_bnd(1)+1e-12],...
    [bbox_size_lower_bnd(2) bbox_size_upper_bnd(2)+1e-12],...
    [bbox_size_lower_bnd(3) bbox_size_upper_bnd(3)+1e-12],0);

% assign constitutive parameters
[idx,epsilon_r,sigma_e,grid_intcon] = intcon_constparams(r,Res,Cnt,Dims,Orients,er,se,0);

% -------------------------------------------------------------------------
%                 Define EM Vars/Constants and Domain Parameters
% -------------------------------------------------------------------------

% Constitutive parameters
Mr = epsilon_r - 1j*sigma_e/(eo*omega); % permittivity
Mc = Mr - 1.0; % susceptibility
OneoverMc = 1.0 ./ Mc; % one over susceptibility

% Domain and unknown parameters
[L,M,N,~] = size(r); % domain size
nD = L*M*N; % number of variables in the system
dx = Res; % voxel size
idxS = find(abs(Mc(:)) > 1e-12); % the indices of the non-air voxel positions
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions for 3 Cartesian components
idxS5 = [idxS; nD+idxS; 2*nD+idxS; 3*nD+idxS; 4*nD+idxS]; % for currents


% -------------------------------------------------------------------------
%                 Specify Excitation and Ground Points
% -------------------------------------------------------------------------

% finding excitation and ground nodes in the spacing between arms

% For first coil

pnt_lft_x=(((n_turns*grd_spc)+width_arm))/2-spc_dg/2;
pnt_rght_x=(((n_turns*grd_spc)+width_arm))/2+spc_dg/2;
dum=1;
pnt_lft=cell(n_arr_elems^2,1);
pnt_rght=cell(n_arr_elems^2,1);
pnt_lft{1}=zeros(round(width_arm/Res)*round(height_arm/Res),3);
pnt_rght{1}=zeros(round(width_arm/Res)*round(height_arm/Res),3);
for kk=1:round(width_arm/Res)
    for ll=1:round(height_arm/Res)
        pnt_lft{1}(dum,1:3)=[pnt_lft_x (2*kk-1)*(0.5*Res) (-len_feed)+(2*ll-1)*(0.5*Res)]; % points on which excitation defined
        pnt_rght{1}(dum,1:3)=[pnt_rght_x (2*kk-1)*(0.5*Res) (-len_feed)+(2*ll-1)*(0.5*Res)]; % points on which ground defined
        dum=dum+1;
    end
end

% For all remaining coils

for kk=2:n_arr_elems^2
    pnt_lft{kk}(:,1)=pnt_lft{1}(:,1) +  shft_vects(kk,1);
    pnt_lft{kk}(:,2)=pnt_lft{1}(:,2) +  shft_vects(kk,2);
    pnt_lft{kk}(:,3)=pnt_lft{1}(:,3) +  shft_vects(kk,3);
    
    pnt_rght{kk}(:,1)=pnt_rght{1}(:,1) + shft_vects(kk,1);
    pnt_rght{kk}(:,2)=pnt_rght{1}(:,2) + shft_vects(kk,2);
    pnt_rght{kk}(:,3)=pnt_rght{1}(:,3) + shft_vects(kk,3);
    
end

% finding the node/current ids of ground and excitation points/unknowns

tend = toc(tini);
disp(['Time for geometry processing ::: ' ,num2str(tend)]);

% temporarily puting points in arrays - to get their node/current ids
num_pnt_tmp=round(width_arm/Res)*round(height_arm/Res);
pnt_exc_dum=zeros(num_pnt_tmp*(n_arr_elems^2),3); % left point
pnt_grnd_dum=zeros(num_pnt_tmp*(n_arr_elems^2),3); % right points

for kk=1:n_arr_elems^2
    pnt_exc_dum((kk-1)*num_pnt_tmp+1:kk*num_pnt_tmp,1:3) = pnt_lft{kk}(:,1:3);
    pnt_grnd_dum((kk-1)*num_pnt_tmp+1:kk*num_pnt_tmp,1:3) = pnt_rght{kk}(:,1:3);
end




