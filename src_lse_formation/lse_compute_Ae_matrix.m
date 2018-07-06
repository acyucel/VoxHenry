function [Ae,nodeid_lft,nodeid_rght,nodeid_wlcond,Ae_only_leaving,Ae_only_entering_bndry] = lse_compute_Ae_matrix(idxS, grid_intcon, L, M, N, dx, pnt_lft, pnt_rght,pnt_well_cond,fl_check_ports)

% Additional codelet for seperate running
%clc; close all; clear all; L=3; M=2; N=2; dx = 0.001;
%grid_tmp = ones(L,M,N); idxS = find(abs(grid_tmp(:)) > 1e-12); clear grid_tmp;
%[grid_intcon] = generategridfrombbox(dx,[0 L*dx],[0 M*dx],[0 N*dx],1);
%fl_meth_comp=2 % methods: 0-> slow , 1-> fast with sort, 2-> fast w/ graph, 3-> fast w/ inspection
%pnt_exc=[]; pnt_grnd=[];

nodeid_4_injectcurr=[];
nodeid_4_grnd=[];
nodeid_4_well_cond=[];

tstart = tic;
fl_profile = 0;
fl_meth_comp = 3;

% constants
num_nonair_cube=length(idxS);

% Inputs for visualization
fl_vis_only_voxels_with_ids=0;
fl_vis_panel_locs_ids_x=0;
fl_vis_panel_locs_ids_y=0;
fl_vis_panel_locs_ids_z=0;
fl_vis_exc_grnd = fl_check_ports;

num_elem_allowed = 5000;
freq_vis = 1;
if (num_nonair_cube > num_elem_allowed)% visualization of boxes is slow when nD > num_elem_allowed
    freq_vis = 25; % every xx element unknown id will be written.
end

% putting port nodes in temporary arrays

% find number of ports
num_ports=size(pnt_lft,1);
if (size(pnt_rght,1) ~= size(pnt_lft,1))
    error('Port definition is wrong !!! Number of faces should 2 for each port')
end

% find number of nodes on each face of each port
num_node_each_port=zeros(num_ports,2);
for kk=1:num_ports
    num_node_each_port(kk,1)=size(pnt_rght{kk},1);
    num_node_each_port(kk,2)=size(pnt_lft{kk},1);
end

% temporarily define the right nodes as ground 
pnt_grnd=zeros(sum(num_node_each_port(:,1)),3); % right points
% temporarily define the left nodes as excitation
pnt_exc=zeros(sum(num_node_each_port(:,2)),3); % left point

% put the nodes in temporary arrays
dum_exc=0;
dum_grnd=0;
for kk=1:num_ports
    pnt_exc(dum_exc+1:dum_exc+num_node_each_port(kk,2),1:3) = pnt_lft{kk}(:,1:3);
    dum_exc=dum_exc+num_node_each_port(kk,2);
    
    pnt_grnd(dum_grnd+1:dum_grnd+num_node_each_port(kk,1),1:3) = pnt_rght{kk}(:,1:3);
    dum_grnd=dum_grnd+num_node_each_port(kk,1);
end

if (fl_meth_comp == 3) % inspection based fast Ae generation
    
    disp('-----------------------------------------------------')
    disp('Generating Ae matrix by inspection...')

    % 1) Get the 'boolean' matrix. Zero if air voxel, non-empty voxel 
    %    sequential number if non-empty. Sequence is x, y, z as per natural
    %    MatLab multi-dimensional matrix index ordering
    tic
    boolean_tens=zeros(L,M,N);
    boolean_tens(idxS)=[1:num_nonair_cube];

    if(fl_profile == 1); disp(['Time for obtaining numbering matrix ::: ', num2str(toc)]); end
    
   
    % 7) Find the panels of each voxel. Voxels are always scan first along x, then y, then z.
    %    However, x panels are numbered the same way, but y panels are ordered y, z, x
    %    and z panels are ordered z, x, y
    
    tic
    
    % Numbering panels along x direction
    %
    % This is easy, as the panels along the x direction follow the same x,y,z voxel order
    % of the voxel scan
    
    % counter for all panels along x
    ind_panel = 1;
    % non-empty voxel counter
    vox_num = 1;
    voxel2panel_x=zeros(num_nonair_cube,2); % 1 left, 2 right panel
    for nn=1:N
        for mm=1:M
            for ll=1:L
                % if non-empty voxel 
                if (boolean_tens(ll,mm,nn) ~=0)
                    % check if next voxel along x (as we are scanning along x)
                    % is empty or not
                    ll_next_one = ll + 1;
                    % if this was the last voxel in the row, there is no next one
                    if (ll_next_one > L) 
                        ll_next_one = -1;
                    % now we know it is not the last one, so we can safely access
                    % the next one and test it for emptyness
                    elseif (boolean_tens(ll_next_one, mm, nn) ==0)
                        ll_next_one = -1;
                    end

                    % assign panels to voxel
                    voxel2panel_x(vox_num,1)=ind_panel; % left panel
                    ind_panel=ind_panel+1;
                    voxel2panel_x(vox_num,2)=ind_panel; % right panel
                    
                    vox_num = vox_num + 1;
                    
                    % if next voxel is empty, or end of the row, increment the 'ind_panel';
                    % otherwise not, as the voxels share the same panel at the interface.
                    if (ll_next_one == -1)
                        ind_panel = ind_panel + 1;
                    end
                end
            end
        end
    end
    num_panel_x=max(max(voxel2panel_x));

    % Numbering panels along y direction
    %
    % This is twisted, as the panels along the y direction follow the y,z,x voxel order
    % which is different from the natural x,y,z order
    
    % counter for all panels along y
    ind_panel = 1;
    voxel2panel_y=zeros(num_nonair_cube,2); % 1 front panel, 2 back panel
    for ll=1:L
        for nn=1:N
          for mm=1:M  
                % if non-empty voxel 
                if (boolean_tens(ll,mm,nn) ~=0)
                    % check if next voxel along y (as we are scanning along y)
                    % is empty or not
                    mm_next_one = mm + 1;
                    % if this was the last voxel in the row, there is no next one
                    if (mm_next_one > M) 
                        mm_next_one = -1;
                    % now we know it is not the last one, so we can safely access
                    % the next one and test it for emptyness
                    elseif (boolean_tens(ll, mm_next_one, nn) ==0)
                        mm_next_one = -1;
                    end

                    % assign panels to voxel
                    %
                    % retrieve the voxel number from the boolean_tens matrix, that has
                    % stored the non-empty voxel sequence number x,y,z 
                    vox_num = boolean_tens(ll,mm,nn);
                    % and use it for assigning the right panel indexes
                    voxel2panel_y(vox_num,1)=ind_panel; % left panel
                    ind_panel=ind_panel+1;
                    voxel2panel_y(vox_num,2)=ind_panel; % right panel
                    
                    % if next voxel is empty, or end of the row, increment the 'ind_panel';
                    % otherwise not, as the voxels share the same panel at the interface.
                    if (mm_next_one == -1)
                        ind_panel = ind_panel + 1;
                    end
                end
            end
        end
    end
    num_panel_y=max(max(voxel2panel_y));
  
    % This is twisted, as the panels along the z direction follow the z,x,y voxel order
    % which is different from the natural x,y,z order
    
    % counter for all panels along z
    ind_panel = 1;
    voxel2panel_z=zeros(num_nonair_cube,2); % 1 front panel, 2 back panel
    for mm=1:M  
        for ll=1:L
            for nn=1:N
                % if non-empty voxel 
                if (boolean_tens(ll,mm,nn) ~=0)
                    % check if next voxel along z (as we are scanning along z)
                    % is empty or not
                    nn_next_one = nn + 1;
                    % if this was the last voxel in the row, there is no next one
                    if (nn_next_one > N) 
                        nn_next_one = -1;
                    % now we know it is not the last one, so we can safely access
                    % the next one and test it for emptyness
                    elseif (boolean_tens(ll, mm, nn_next_one) ==0)
                        nn_next_one = -1;
                    end

                    % assign panels to voxel
                    %
                    % retrieve the voxel number from the boolean_tens matrix, that has
                    % stored the non-empty voxel sequence number x,y,z 
                    vox_num = boolean_tens(ll,mm,nn);
                    % and use it for assigning the right panel indexes
                    voxel2panel_z(vox_num,1)=ind_panel; % left panel
                    ind_panel=ind_panel+1;
                    voxel2panel_z(vox_num,2)=ind_panel; % right panel
                    
                    % if next voxel is empty, or end of the row, increment the 'ind_panel';
                    % otherwise not, as the voxels share the same panel at the interface.
                    if (nn_next_one == -1)
                        ind_panel = ind_panel + 1;
                    end
                end
            end
        end
    end
    num_panel_z=max(max(voxel2panel_z));
  
    %
    % here we continue as with the 'old' Graph-based routine
    %    
    
    voxel2panel_y = num_panel_x + voxel2panel_y;
    
    voxel2panel_z = (num_panel_x+num_panel_y) + voxel2panel_z;
    
    all_panels_ids=[voxel2panel_x; voxel2panel_y; voxel2panel_z];
    
    if(fl_profile == 1); disp(['Time for indexing panels ::: ', num2str(toc)]); end;
    
    % 8) forming Ae matrix
    
    tic
    num_nodes=num_panel_x+num_panel_y+num_panel_z;
    
    % enter +1, leaves -1; first entry for leaving, second entry for entering
    const_lin=1/2;
    %const_lin=dx/2;
    
    sp_mat_inds=zeros(16*num_nonair_cube,3);
    sp_mat_inds_only_leaving_currs=zeros(8*num_nonair_cube,3);
    for kk=1:3*num_nonair_cube % pertinent to Jx, Jy, and Jz
        sp_mat_inds(kk,1:3)= [all_panels_ids(kk,1) kk -1];
        sp_mat_inds(3*num_nonair_cube+kk,1:3)= [all_panels_ids(kk,2) kk 1];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(kk,1:3)= [all_panels_ids(kk,1) kk -1];
    end
    
    dum=1;
    for kk=3*num_nonair_cube+1:4*num_nonair_cube % pertinent to Jfx
        sp_mat_inds(6*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
        sp_mat_inds(7*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,2) kk const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(3*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
        dum=dum+1;
    end
    
    dum=1;
    for kk=3*num_nonair_cube+1:4*num_nonair_cube % pertinent to Jfy
        sp_mat_inds(8*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk -const_lin];
        sp_mat_inds(9*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,2) kk -const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(4*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk -const_lin];
        dum=dum+1;
    end
    
    % for space diagonal currents
    dum=1;
    for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsx
        sp_mat_inds(10*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
        sp_mat_inds(11*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,2) kk const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(5*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
        dum=dum+1;
    end
    
    dum=1;
    for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsy
        sp_mat_inds(12*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk const_lin];
        sp_mat_inds(13*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,2) kk const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(6*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk const_lin];
        dum=dum+1;
    end
    
    dum=1;
    for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsz
        sp_mat_inds(14*num_nonair_cube+dum,1:3)= [all_panels_ids(2*num_nonair_cube+dum,1) kk -2*const_lin];
        sp_mat_inds(15*num_nonair_cube+dum,1:3)= [all_panels_ids(2*num_nonair_cube+dum,2) kk -2*const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(7*num_nonair_cube+dum,1:3)= [all_panels_ids(2*num_nonair_cube+dum,1) kk -2*const_lin];
        dum=dum+1;
    end
    
    %Ae=sparse(sp_mat_inds(:,1),sp_mat_inds(:,2),sp_mat_inds(:,3));
    Ae=sparse(sp_mat_inds(:,1),sp_mat_inds(:,2),sp_mat_inds(:,3),num_nodes,5*num_nonair_cube);
    
    infomem1 = whos('Ae');
    memestimated = (infomem1.bytes)/(1024*1024);
    disp(['Memory for Ae matrix (MB)::' , num2str(memestimated)]);
    
    % Additional data structure for visualization of currents
    
    % We need two data structures:
    % (1) Ae matrix for leaving currents
    Ae_only_leaving=sparse(sp_mat_inds_only_leaving_currs(:,1),sp_mat_inds_only_leaving_currs(:,2),sp_mat_inds_only_leaving_currs(:,3),num_nodes,5*num_nonair_cube);
    
    % (2) Ae matrix for entering currents to boundary nodes
    % Finding the nodes to which current enters
    dum_vect=zeros(size(Ae,2),1);
    dum_vect(1:size(Ae,2)/5*3)=1;
    dum_res=Ae*dum_vect; % ids of boundary nodes -1:exiting, +1:entering
    
    nodes_only_curr_enter=find(dum_res == 1);
    sp_mat_dum=Ae(nodes_only_curr_enter,:);
    [rows,cols,vals] = find(sp_mat_dum);
    for kk=1:length(rows) %correct rows
        rows(kk)=nodes_only_curr_enter(rows(kk));
    end
    
    Ae_only_entering_bndry=sparse(rows,cols,vals,num_nodes,5*num_nonair_cube);
    
    if(fl_profile == 1); disp(['Time for generating sparse Ae matrix ::: ', num2str(toc)]);end;
    
    % The following is for visualization!!!
    % retrieving the centers of non-air voxels
    tic
    dum=1;
    xy_curr_ids_locs=zeros(3*num_nonair_cube,3);
    for mm=1:N
        for ll=1:M
            for kk=1:L
                if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                        grid_intcon(kk,ll,mm,3) ~=0 ) % This is dangerous!
                    coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                    xy_curr_ids_locs(dum,1:3)=coor_tmp;
                    xy_curr_ids_locs(num_nonair_cube+dum,1:3)=coor_tmp;
                    xy_curr_ids_locs(2*num_nonair_cube+dum,1:3)=coor_tmp;
                    dum=dum+1;
                end
            end
        end
    end
    if(fl_profile == 1); disp(['Time for retrieving centers of non-air voxels ::: ', num2str(toc)]);end;
    
    
    % finding locations of panels enclosing x and y currents
    tic
    all_panel_locs=zeros(3*num_nonair_cube,6);
    for kk=1:num_nonair_cube
        %Jx currents
        all_panel_locs(kk,1:3)=[xy_curr_ids_locs(kk,1)-dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
        all_panel_locs(kk,4:6)=[xy_curr_ids_locs(kk,1)+dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
        %Jy currents
        ind=num_nonair_cube+kk;
        all_panel_locs(ind,1:3)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)-dx*0.5 xy_curr_ids_locs(ind,3)];
        all_panel_locs(ind,4:6)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)+dx*0.5 xy_curr_ids_locs(ind,3)];
        %Jz currents
        ind2=2*num_nonair_cube+kk;
        all_panel_locs(ind2,1:3)=[xy_curr_ids_locs(ind2,1) xy_curr_ids_locs(ind2,2) xy_curr_ids_locs(ind2,3)-dx*0.5];
        all_panel_locs(ind2,4:6)=[xy_curr_ids_locs(ind2,1) xy_curr_ids_locs(ind2,2) xy_curr_ids_locs(ind2,3)+dx*0.5];
    end
    if(fl_profile == 1); disp(['Time for finding locations of surface panels ::: ', num2str(toc)]); end
    

    
    fl_vis_numbering=0;
    if (fl_vis_numbering == 1)
        figure;
        set(gca,'FontSize',24);
        xd = grid_intcon(:,:,:,1);yd = grid_intcon(:,:,:,2);zd = grid_intcon(:,:,:,3);
        h=plot3(xd(:), yd(:), zd(:), 'r*');
        set(h,'MarkerSize',10); xlabel('x'); ylabel('y'); zlabel('z');set(gca,'FontSize',24);
        hold on
        for mm=1:N
            for ll=1:M
                for kk=1:L
                    coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                    %h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(numbering_x_y_z(kk,ll,mm,1))); % x numbering
                    h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(numbering_x_y_z(kk,ll,mm,2))); % y numbering
                    %h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(numbering_x_y_z(kk,ll,mm,3))); % z numbering
                    set(h,'FontSize',24)
                    
                    
                end
            end
        end
        xlabel('x');ylabel('y'); set(gca,'FontSize',24);
        grid on
        axis tight
        grid on
        view(2)
    end

    disp('Done... Generating Ae matrix by inspection...')
    disp('-----------------------------------------------------')
    
elseif (fl_meth_comp == 2) % graph based fast Ae generation
    
    disp('-----------------------------------------------------')
    disp('Generating Ae matrix with graphs...')
    disp('This routine uses graph toolbox of matlab')
    
    % 1) Get the boolean tensor
    tic
    boolean_tens=zeros(L,M,N);
    boolean_tens(idxS)=1;
    if(fl_profile == 1); disp(['Time for obtaining boolean matrix ::: ', num2str(toc)]); end
    
    % 2) Number the non-empty voxels (or nodes in graph) along x, y, and z
    % directions. Three different number is used to ensure the correct ordering
    % of panels along x, y, and z directions in the graphs.
    
    tic
    % 2a) ijk to ind tensor and ind to ijk tensor for numbering along x(x_num_vect)
    unkids_ijk_to_ind=zeros(L,M,N);
    x_num_vect=zeros(num_nonair_cube,3);
    
    dum=0;
    for mm=1:N
        for ll=1:M
            for kk=1:L
                if (boolean_tens(kk,ll,mm) ~=0)
                    dum=dum+1;
                    unkids_ijk_to_ind(kk,ll,mm)=dum;
                    x_num_vect(dum,1:3)=[kk ll mm];
                end
            end
        end
    end
    
    %2b) ind to ijk tensor for numbering along y(y_num_vect)
    y_num_vect=zeros(num_nonair_cube,3);
    dum=0;
    for kk=1:L
        for mm=1:N
            for ll=1:M
                if (boolean_tens(kk,ll,mm) ~=0)
                    dum=dum+1;
                    y_num_vect(dum,1:3)=[kk ll mm];
                end
            end
        end
    end
    
    %2c) ind to ijk tensor for numbering along z(z_num_vect)
    z_num_vect=zeros(num_nonair_cube,3);
    dum=0;
    for ll=1:M
        for kk=1:L
            for mm=1:N
                if (boolean_tens(kk,ll,mm) ~=0)
                    dum=dum+1;
                    z_num_vect(dum,1:3)=[kk ll mm];
                end
            end
        end
    end
    
    % 2d) ijk tensor storing the numbering along x, y, and z directions
    numbering_x_y_z=zeros(L,M,N,3);
    for kk=1:size(x_num_vect,1)
        numbering_x_y_z(x_num_vect(kk,1),x_num_vect(kk,2),x_num_vect(kk,3),1)=kk; % x numbering
        numbering_x_y_z(y_num_vect(kk,1),y_num_vect(kk,2),y_num_vect(kk,3),2)=kk; % y numbering
        numbering_x_y_z(z_num_vect(kk,1),z_num_vect(kk,2),z_num_vect(kk,3),3)=kk; % z numbering
    end
    
    %2e) mapping from y and z numbering to x numbering
    
    tic
    dum=0;
    map_from_y2x_numbering=zeros(num_nonair_cube,1);
    map_from_z2x_numbering=zeros(num_nonair_cube,1);
    %for kk=1:L
    %    for mm=1:N
    %        for ll=1:M
    %            if (boolean_tens(kk,ll,mm) ~=0)
    %                dum=dum+1;
    %                map_from_y2x_numbering(dum)=numbering_x_y_z(y_num_vect(dum,1),y_num_vect(dum,2),y_num_vect(dum,3),1);
    %                map_from_z2x_numbering(dum)=numbering_x_y_z(z_num_vect(dum,1),z_num_vect(dum,2),z_num_vect(dum,3),1);
    %             end
    %        end
    %    end
    %end
    % faster code 
    for dum=1:num_nonair_cube
        map_from_y2x_numbering(dum)=numbering_x_y_z(y_num_vect(dum,1),y_num_vect(dum,2),y_num_vect(dum,3),1);
        map_from_z2x_numbering(dum)=numbering_x_y_z(z_num_vect(dum,1),z_num_vect(dum,2),z_num_vect(dum,3),1);
    end
   
    if(fl_profile == 1); disp(['Time for getting numbering matrices ::: ', num2str(toc)]);end
    clear y_num_vect z_num_vect
    
    % 3) Finding the near neighbors of each non-empty voxel
    tic
    temp_crit=1; % box_diff
    tola=1e-12;
    nn_nums_tens=zeros(L,M,N,3); % number of near neighbors of each voxel
    nn_x_neigh_ind_tens=zeros(num_nonair_cube,2); % 1-> left, 2-> right neighbors
    nn_y_neigh_ind_tens=zeros(num_nonair_cube,2); % 1-> front, 2-> back neighbors
    nn_z_neigh_ind_tens=zeros(num_nonair_cube,2); % 1-> bottom, 2-> up neighbors
    
    for mm=1:N%1 % z variation
        for ll=1:M % y variation
            for kk=1:L % x variation
                if (boolean_tens(kk,ll,mm) ~=0) %non-air voxel (source)
                    
                    lowbnd_indi=kk-temp_crit;
                    upbnd_indi=kk+temp_crit;
                    
                    lowbnd_indj=ll-temp_crit;
                    upbnd_indj=ll+temp_crit;
                    
                    lowbnd_indk=mm-temp_crit;
                    upbnd_indk=mm+temp_crit;
                    
                    dum=0; % counter for near neighbors of each voxel with index (kk,ll,mm)
                    % counters for near neighbor along x,y,and z directions
                    dum_nn_x_cnt=0; dum_nn_y_cnt=0; dum_nn_z_cnt=0;
                    for cc=lowbnd_indk:upbnd_indk % search along z direction
                        
                        if (cc >= 1 && cc <= N)
                            
                            for bb=lowbnd_indj:upbnd_indj % search along y direction
                                
                                if (bb >= 1 && bb <= M)
                                    
                                    for aa=lowbnd_indi:upbnd_indi % search along x direction
                                        
                                        if (aa >= 1 && aa<= L )
                                            
                                            if (boolean_tens(aa,bb,cc)==1) %non-air voxel (observer)
                                                
                                                diffind=abs(aa-kk)+abs(bb-ll)+abs(cc-mm);
                                                
                                                if( diffind < 2-tola && diffind > 0+tola)% corner and self neighbors
                                                    
                                                    if (abs(abs(aa-kk)-1) < tola) % x neighbor
                                                        if (sign(aa-kk) < 0) % left neighbor
                                                            nn_x_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,1),1) = numbering_x_y_z(aa,bb,cc,1);
                                                        else % right neighbor
                                                            nn_x_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,1),2) = numbering_x_y_z(aa,bb,cc,1);
                                                        end
                                                        dum_nn_x_cnt=dum_nn_x_cnt+1;
                                                    elseif (abs(abs(bb-ll)-1) < tola) % y neighbor
                                                        if (sign(bb-ll) < 0) % front neighbor
                                                            nn_y_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,2),1) = numbering_x_y_z(aa,bb,cc,2);
                                                        else % back neighbor
                                                            nn_y_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,2),2) = numbering_x_y_z(aa,bb,cc,2);
                                                        end
                                                        dum_nn_y_cnt=dum_nn_y_cnt+1;
                                                    elseif (abs(abs(cc-mm)-1) < tola) % z neighbor
                                                        if (sign(cc-mm) < 0) % bottom neighbor
                                                            nn_z_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,3),1) = numbering_x_y_z(aa,bb,cc,3);
                                                        else % up neighbor
                                                            nn_z_neigh_ind_tens(numbering_x_y_z(kk,ll,mm,3),2) = numbering_x_y_z(aa,bb,cc,3);
                                                        end
                                                        dum_nn_z_cnt=dum_nn_z_cnt+1;
                                                    end
                                                    dum=dum+1;
                                                    
                                                end
                                                
                                            end
                                            
                                        end
                                        
                                    end
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                    nn_nums_tens(kk,ll,mm,1)=dum_nn_x_cnt;
                    nn_nums_tens(kk,ll,mm,2)=dum_nn_y_cnt;
                    nn_nums_tens(kk,ll,mm,3)=dum_nn_z_cnt;
                    
                end
                
            end
            
        end
        
    end
    if(fl_profile == 1); disp(['Time for determining neighbors of voxels ::: ', num2str(toc)]); end
    
    % 4) Create adjacency matrix for the nodes with numbering along x, y,and
    % z directions
    tic
    
    % Note: if we create the indices for sparse adjacency matrix and then
    % allocate the sparse matrix with sparse command, we see that matrix
    % doesn't have the dimension of (num_nonair_cubexnum_nonair_cube)
    % (if the last elements are zeros) We should specify the dimensions of
    % sparse matrices in sparse command
    
    % adjacency matrix for x directed panels
    inds_dum=zeros(sum(sum(sum(nn_nums_tens(:,:,:,1)))),3);
    dum=0;
    for kk=1:size(nn_x_neigh_ind_tens,1)
        if(nn_x_neigh_ind_tens(kk,1) > 0)
            dum=dum+1;
            inds_dum(dum,1:3)=[kk nn_x_neigh_ind_tens(kk,1) 1];
        end
        if(nn_x_neigh_ind_tens(kk,2) > 0)
            dum=dum+1;
            inds_dum(dum,1:3)=[kk nn_x_neigh_ind_tens(kk,2) 1];
        end
    end
    
    adj_mat_x=sparse(inds_dum(:,1),inds_dum(:,2),inds_dum(:,3),num_nonair_cube,num_nonair_cube);
    
    if(isempty(adj_mat_x) == 1)
        adj_mat_x=sparse(num_nonair_cube,num_nonair_cube);
    end
    
    % adjacency matrix for y directed panels
    
    inds_dum=zeros(sum(sum(sum(nn_nums_tens(:,:,:,2)))),3);
    dum=0;
    for kk=1:size(nn_y_neigh_ind_tens,1)
        if(nn_y_neigh_ind_tens(kk,1) > 0)
            dum=dum+1;
            inds_dum(dum,1:3)=[kk nn_y_neigh_ind_tens(kk,1) 1];
        end
        if(nn_y_neigh_ind_tens(kk,2) > 0)
            dum=dum+1;
            inds_dum(dum,1:3)=[kk nn_y_neigh_ind_tens(kk,2) 1];
        end
    end
    
    adj_mat_y=sparse(inds_dum(:,1),inds_dum(:,2),inds_dum(:,3),num_nonair_cube,num_nonair_cube);
    
    if(isempty(adj_mat_y) == 1)
        adj_mat_y=sparse(num_nonair_cube,num_nonair_cube);
    end
    
    % adjacency matrix for z directed panels
    
    inds_dum=zeros(sum(sum(sum(nn_nums_tens(:,:,:,3)))),3);
    dum=0;
    for kk=1:size(nn_z_neigh_ind_tens,1)
        if(nn_z_neigh_ind_tens(kk,1) > 0)
            dum=dum+1;
            inds_dum(dum,1:3)=[kk nn_z_neigh_ind_tens(kk,1) 1];
        end
        if(nn_z_neigh_ind_tens(kk,2) > 0)
            dum=dum+1;
            inds_dum(dum,1:3)=[kk nn_z_neigh_ind_tens(kk,2) 1];
        end
    end
    
    adj_mat_z=sparse(inds_dum(:,1),inds_dum(:,2),inds_dum(:,3),num_nonair_cube,num_nonair_cube);
    
    if(isempty(adj_mat_z) == 1)
        adj_mat_z=sparse(num_nonair_cube,num_nonair_cube);
    end
    
    if(fl_profile == 1); disp(['Time for getting adjacency matrices ::: ', num2str(toc)]); end;
    
    %disp('Are adjacency matrices for x,y,and z symmetric?')
    %[issymmetric(adj_mat_x) issymmetric(adj_mat_y) issymmetric(adj_mat_z)]
    
    % 5) Create graphs using adjacency matrices
    tic
    G_x = graph(adj_mat_x);
    G_y = graph(adj_mat_y);
    G_z = graph(adj_mat_z);
    if(fl_profile == 1); disp(['Time for creating graphs ::: ', num2str(toc)]); end;
    
    % 6) Find the connected elements in each graph
    tic
    bins_x=conncomp(G_x,'OutputForm','cell');
    bins_y=conncomp(G_y,'OutputForm','cell');
    bins_z=conncomp(G_z,'OutputForm','cell');
    if(fl_profile == 1); disp(['Time for finding conn elements in graphs ::: ', num2str(toc)]);end;
    
    % 7) Find the panels of each voxel (due to x numbering)
    tic
    % Numbering panels along x direction
    ind_panel=0; % counter for all panels along x
    voxel2panel_x=zeros(num_nonair_cube,2); % 1 left, 2 right panel
    for kk=1:length(bins_x)
        
        for ll=1:length(bins_x{kk})
            
            if (ll == 1) % beginning of connected elements, put a panel
                ind_panel=ind_panel+1;
            end
            
            unk_id=bins_x{kk}(ll);
            
            voxel2panel_x(unk_id,1)=ind_panel; % left panel
            
            ind_panel=ind_panel+1; % this will automatically put panel at the end of connected elements
            
            voxel2panel_x(unk_id,2)=ind_panel; % right panel
        end
        
    end
    num_panel_x=max(max(voxel2panel_x));
    
    % Numbering panels along y direction
    ind_panel=0; % counter for all panels along y
    voxel2panel_y=zeros(num_nonair_cube,2); % 1 front, 2 back
    for kk=1:length(bins_y)
        
        for ll=1:length(bins_y{kk})
            
            if (ll == 1) % beginning of connected elements, put a panel
                ind_panel=ind_panel+1;
            end
            
            unk_id=map_from_y2x_numbering(bins_y{kk}(ll));
            
            voxel2panel_y(unk_id,1)=ind_panel; % front panel
            
            ind_panel=ind_panel+1; % this will automatically put panel at the end of connected elements
            
            voxel2panel_y(unk_id,2)=ind_panel; % back panel
        end
        
    end
    num_panel_y=max(max(voxel2panel_y));
    
    % Numbering panels along z direction
    ind_panel=0; % counter for all panels along z
    voxel2panel_z=zeros(num_nonair_cube,2); % 1 bottom, 2 up
    for kk=1:length(bins_z)
        
        for ll=1:length(bins_z{kk})
            
            if (ll == 1) % beginning of connected elements, put a panel
                ind_panel=ind_panel+1;
            end
            
            unk_id=map_from_z2x_numbering(bins_z{kk}(ll));
            
            voxel2panel_z(unk_id,1)=ind_panel; % front panel
            
            ind_panel=ind_panel+1; % this will automatically put panel at the end of connected elements
            
            voxel2panel_z(unk_id,2)=ind_panel; % back panel
        end
        
    end
    num_panel_z=max(max(voxel2panel_z));
    
    voxel2panel_y = num_panel_x + voxel2panel_y;
    
    voxel2panel_z = (num_panel_x+num_panel_y) + voxel2panel_z;
    
    all_panels_ids=[voxel2panel_x; voxel2panel_y; voxel2panel_z];
    if(fl_profile == 1); disp(['Time for indexing panels ::: ', num2str(toc)]); end;
    
    % 8) forming Ae matrix
    
    tic
    num_nodes=num_panel_x+num_panel_y+num_panel_z;
    
    % enter +1, leaves -1; first entry for leaving, second entry for entering
    const_lin=1/2;
    %const_lin=dx/2;
    
    sp_mat_inds=zeros(16*num_nonair_cube,3);
    sp_mat_inds_only_leaving_currs=zeros(8*num_nonair_cube,3);
    for kk=1:3*num_nonair_cube % pertinent to Jx, Jy, and Jz
        sp_mat_inds(kk,1:3)= [all_panels_ids(kk,1) kk -1];
        sp_mat_inds(3*num_nonair_cube+kk,1:3)= [all_panels_ids(kk,2) kk 1];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(kk,1:3)= [all_panels_ids(kk,1) kk -1];
    end
    
    dum=1;
    for kk=3*num_nonair_cube+1:4*num_nonair_cube % pertinent to Jfx
        sp_mat_inds(6*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
        sp_mat_inds(7*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,2) kk const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(3*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
        dum=dum+1;
    end
    
    dum=1;
    for kk=3*num_nonair_cube+1:4*num_nonair_cube % pertinent to Jfy
        sp_mat_inds(8*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk -const_lin];
        sp_mat_inds(9*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,2) kk -const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(4*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk -const_lin];
        dum=dum+1;
    end
    
    % for space diagonal currents
    dum=1;
    for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsx
        sp_mat_inds(10*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
        sp_mat_inds(11*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,2) kk const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(5*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
        dum=dum+1;
    end
    
    dum=1;
    for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsy
        sp_mat_inds(12*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk const_lin];
        sp_mat_inds(13*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,2) kk const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(6*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk const_lin];
        dum=dum+1;
    end
    
    dum=1;
    for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsz
        sp_mat_inds(14*num_nonair_cube+dum,1:3)= [all_panels_ids(2*num_nonair_cube+dum,1) kk -2*const_lin];
        sp_mat_inds(15*num_nonair_cube+dum,1:3)= [all_panels_ids(2*num_nonair_cube+dum,2) kk -2*const_lin];
        % the following is added for current vis.
        sp_mat_inds_only_leaving_currs(7*num_nonair_cube+dum,1:3)= [all_panels_ids(2*num_nonair_cube+dum,1) kk -2*const_lin];
        dum=dum+1;
    end
    
    %Ae=sparse(sp_mat_inds(:,1),sp_mat_inds(:,2),sp_mat_inds(:,3));
    Ae=sparse(sp_mat_inds(:,1),sp_mat_inds(:,2),sp_mat_inds(:,3),num_nodes,5*num_nonair_cube);
    
    infomem1 = whos('Ae');
    memestimated = (infomem1.bytes)/(1024*1024);
    disp(['Memory for Ae matrix (MB)::' , num2str(memestimated)]);
    
    % Additional data structure for visualization of currents
    
    % We need two data structures:
    % (1) Ae matrix for leaving currents
    Ae_only_leaving=sparse(sp_mat_inds_only_leaving_currs(:,1),sp_mat_inds_only_leaving_currs(:,2),sp_mat_inds_only_leaving_currs(:,3),num_nodes,5*num_nonair_cube);
    
    % (2) Ae matrix for entering currents to boundary nodes
    % Finding the nodes to which current enters
    dum_vect=zeros(size(Ae,2),1);
    dum_vect(1:size(Ae,2)/5*3)=1;
    dum_res=Ae*dum_vect; % ids of boundary nodes -1:exiting, +1:entering
    
    nodes_only_curr_enter=find(dum_res == 1);
    sp_mat_dum=Ae(nodes_only_curr_enter,:);
    [rows,cols,vals] = find(sp_mat_dum);
    for kk=1:length(rows) %correct rows
        rows(kk)=nodes_only_curr_enter(rows(kk));
    end
    
    Ae_only_entering_bndry=sparse(rows,cols,vals,num_nodes,5*num_nonair_cube);
    
    if(fl_profile == 1); disp(['Time for generating sparse Ae matrix ::: ', num2str(toc)]);end;
    
    % The following is for visualization!!!
    % retrieving the centers of non-air voxels
    tic
    dum=1;
    xy_curr_ids_locs=zeros(3*num_nonair_cube,3);
    for mm=1:N
        for ll=1:M
            for kk=1:L
                if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                        grid_intcon(kk,ll,mm,3) ~=0 ) % This is dangerous!
                    coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                    xy_curr_ids_locs(dum,1:3)=coor_tmp;
                    xy_curr_ids_locs(num_nonair_cube+dum,1:3)=coor_tmp;
                    xy_curr_ids_locs(2*num_nonair_cube+dum,1:3)=coor_tmp;
                    dum=dum+1;
                end
            end
        end
    end
    if(fl_profile == 1); disp(['Time for retrieving centers of non-air voxels ::: ', num2str(toc)]);end;
    
    
    % finding locations of panels enclosing x and y currents
    tic
    all_panel_locs=zeros(3*num_nonair_cube,6);
    for kk=1:num_nonair_cube
        %Jx currents
        all_panel_locs(kk,1:3)=[xy_curr_ids_locs(kk,1)-dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
        all_panel_locs(kk,4:6)=[xy_curr_ids_locs(kk,1)+dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
        %Jy currents
        ind=num_nonair_cube+kk;
        all_panel_locs(ind,1:3)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)-dx*0.5 xy_curr_ids_locs(ind,3)];
        all_panel_locs(ind,4:6)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)+dx*0.5 xy_curr_ids_locs(ind,3)];
        %Jz currents
        ind2=2*num_nonair_cube+kk;
        all_panel_locs(ind2,1:3)=[xy_curr_ids_locs(ind2,1) xy_curr_ids_locs(ind2,2) xy_curr_ids_locs(ind2,3)-dx*0.5];
        all_panel_locs(ind2,4:6)=[xy_curr_ids_locs(ind2,1) xy_curr_ids_locs(ind2,2) xy_curr_ids_locs(ind2,3)+dx*0.5];
    end
    if(fl_profile == 1); disp(['Time for finding locations of surface panels ::: ', num2str(toc)]); end
    
    % visualize graphs along x, y, or z directions
    fl_vis_graphs=0;
    if(fl_vis_graphs == 1)
        figure; plot(G_x);title('Graph for panels along x direction')
        figure; plot(G_y);title('Graph for panels along y direction')
        figure; plot(G_z);title('Graph for panels along z direction')
    end
    
    fl_vis_numbering=0;
    if (fl_vis_numbering == 1)
        figure;
        set(gca,'FontSize',24);
        xd = grid_intcon(:,:,:,1);yd = grid_intcon(:,:,:,2);zd = grid_intcon(:,:,:,3);
        h=plot3(xd(:), yd(:), zd(:), 'r*');
        set(h,'MarkerSize',10); xlabel('x'); ylabel('y'); zlabel('z');set(gca,'FontSize',24);
        hold on
        for mm=1:N
            for ll=1:M
                for kk=1:L
                    coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                    %h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(numbering_x_y_z(kk,ll,mm,1))); % x numbering
                    h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(numbering_x_y_z(kk,ll,mm,2))); % y numbering
                    %h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(numbering_x_y_z(kk,ll,mm,3))); % z numbering
                    set(h,'FontSize',24)
                    
                    
                end
            end
        end
        xlabel('x');ylabel('y'); set(gca,'FontSize',24);
        grid on
        axis tight
        grid on
        view(2)
    end

    disp('Done... Generating Ae matrix with graphs...')
    disp('-----------------------------------------------------')

elseif (fl_meth_comp == 1) % sorting based fast Ae generation
    
    disp('-----------------------------------------------------')
    disp('Fast Ae matrix generation with sorting is ON!!! With this routine,')
    disp('correct ordering of nodes is NOT ensured!!! This is because')
    disp('sortrow in uniquetol command of matlab does not always sort the elements')
    disp('in 2nd column of coordinates array (i.e., y-coordinates) when the elements')
    disp('in 1st column of coordinates array (i.e., x-coordinates) are the same.')
    disp('-----------------------------------------------------')
    
    % retrieving the centers of non-air voxels
    tic
    dum=1;
    xy_curr_ids_locs=zeros(2*num_nonair_cube,3);
    for mm=1:N
        for ll=1:M
            for kk=1:L
                if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                        grid_intcon(kk,ll,mm,3) ~=0 ) % This is dangerous!
                    coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                    xy_curr_ids_locs(dum,1:3)=coor_tmp;
                    xy_curr_ids_locs(num_nonair_cube+dum,1:3)=coor_tmp;
                    dum=dum+1;
                end
            end
        end
    end
    if(fl_profile == 1); disp(['Time for retrieving centers of non-air voxels ::: ', num2str(toc)]); end;
    
    
    % finding locations of panels enclosing x and y currents
    tic
    all_panel_locs=zeros(2*num_nonair_cube,6);
    for kk=1:num_nonair_cube
        %Jx currents
        all_panel_locs(kk,1:3)=[xy_curr_ids_locs(kk,1)-dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
        all_panel_locs(kk,4:6)=[xy_curr_ids_locs(kk,1)+dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
        %Jy currents
        ind=num_nonair_cube+kk;
        all_panel_locs(ind,1:3)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)-dx*0.5 xy_curr_ids_locs(ind,3)];
        all_panel_locs(ind,4:6)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)+dx*0.5 xy_curr_ids_locs(ind,3)];
    end
    if(fl_profile == 1); disp(['Time for finding locations of surface panels ::: ', num2str(toc)]); end;
    
    % Finding unique panels in a fast way
    tic
    unique_panel_coors_x=[all_panel_locs(1:num_nonair_cube,1:3);all_panel_locs(1:num_nonair_cube,4:6)];
    unique_panel_coors_y=[all_panel_locs(num_nonair_cube+1:2*num_nonair_cube,1:3);all_panel_locs(num_nonair_cube+1:2*num_nonair_cube,4:6)];
    
    tola=1e-12;
    [unique_panel_coors_x,idx_map] = uniquetol(unique_panel_coors_x,tola,'ByRows',true,'OutputAllIndices',true);
    [unique_panel_coors_y,idy_map] = uniquetol(unique_panel_coors_y,tola,'ByRows',true,'OutputAllIndices',true);
    
    if(fl_profile == 1); disp(['Time for finding unique surface panels ::: ', num2str(toc)]); end;
    
    % Assigning ids of unique panels
    tic
    all_panels_ids_x=zeros(num_nonair_cube,2);
    for kk=1:size(idx_map)
        for ll=1:length(idx_map{kk})
            if (idx_map{kk}(ll) > num_nonair_cube)
                ent_ind=idx_map{kk}(ll)-num_nonair_cube;
            else
                ent_ind=idx_map{kk}(ll);
            end
            if (all_panels_ids_x(ent_ind,1) == 0)
                all_panels_ids_x(ent_ind,1)=kk;
            else
                all_panels_ids_x(ent_ind,2)=kk;
            end
        end
    end
    
    all_panels_ids_y=zeros(num_nonair_cube,2);
    for kk=1:size(idy_map)
        for ll=1:length(idy_map{kk})
            if (idy_map{kk}(ll) > num_nonair_cube)
                ent_ind=idy_map{kk}(ll)-num_nonair_cube;
            else
                ent_ind=idy_map{kk}(ll);
            end
            if (all_panels_ids_y(ent_ind,1) == 0)
                all_panels_ids_y(ent_ind,1)=kk;
            else
                all_panels_ids_y(ent_ind,2)=kk;
            end
        end
    end
    
    max_id_x=max(max(all_panels_ids_x));
    
    all_panels_ids_y = max_id_x + all_panels_ids_y;
    
    all_panels_ids=[all_panels_ids_x;all_panels_ids_y];
    
    if(fl_profile == 1); disp(['Time for finding the ids of unique panels ::: ', num2str(toc)]); end;
    
    % Generate Ae nodal incidence matrix
    tic
    num_node_enc_x=size(unique_panel_coors_x,1); % # of panels enclosing Jx
    num_node_enc_y=size(unique_panel_coors_y,1); % # of panels enclosing Jy
    num_nodes=num_node_enc_x+num_node_enc_y;
    
    %Ae = spalloc(num_nodes,3*num_nonair_cube,8*num_nonair_cube);
    %Ae =zeros(num_nodes,3*num_nonair_cube); % each non air cube has Jx, Jy, Jd
    
    % enter +1, leaves -1; first entry for leaving, second entry for entering
    const_lin=1/2;
    %const_lin=dx/2;
    
    sp_mat_inds=zeros(8*num_nonair_cube,3);
    for kk=1:2*num_nonair_cube
        sp_mat_inds(kk,1:3)= [all_panels_ids(kk,1) kk -1];
        sp_mat_inds(2*num_nonair_cube+kk,1:3)= [all_panels_ids(kk,2) kk 1];
    end
    
    dum=1;
    for kk=2*num_nonair_cube+1:3*num_nonair_cube
        sp_mat_inds(4*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
        sp_mat_inds(5*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,2) kk const_lin];
        dum=dum+1;
    end
    
    dum=1;
    for kk=2*num_nonair_cube+1:3*num_nonair_cube
        sp_mat_inds(6*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk -const_lin];
        sp_mat_inds(7*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,2) kk -const_lin];
        dum=dum+1;
    end
    
    Ae=sparse(sp_mat_inds(:,1),sp_mat_inds(:,2),sp_mat_inds(:,3));
    
    infomem1 = whos('Ae');
    memestimated = (infomem1.bytes)/(1024*1024);
    disp(['Memory for Ae matrix (MB)::' , num2str(memestimated)]);
    
    if(fl_profile == 1); disp(['Time for generating sparse Ae matrix ::: ', num2str(toc)]); end;
    
elseif (fl_meth_comp == 0) % sorting based slow Ae generation
    
    disp('-----------------------------------------------------')
    disp('Slow Ae matrix generation is ON!!! With this routine,')
    disp('correct ordering of nodes is ensured!!!')
    disp('-----------------------------------------------------')
    
    % retrieving the centers of non-air voxels
    tic
    dum=1;
    xy_curr_ids_locs=zeros(2*num_nonair_cube,3);
    for mm=1:N
        for ll=1:M
            for kk=1:L
                if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                        grid_intcon(kk,ll,mm,3) ~=0 ) % This is dangerous!
                    coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                    xy_curr_ids_locs(dum,1:3)=coor_tmp;
                    xy_curr_ids_locs(num_nonair_cube+dum,1:3)=coor_tmp;
                    dum=dum+1;
                end
            end
        end
    end
    if(fl_profile == 1); disp(['Time for retrieving centers of non-air voxels ::: ', num2str(toc)]); end;
    
    % finding locations of panels enclosing x and y currents
    tic
    all_panel_locs=zeros(2*num_nonair_cube,6);
    for kk=1:num_nonair_cube
        %Jx currents
        all_panel_locs(kk,1:3)=[xy_curr_ids_locs(kk,1)-dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
        all_panel_locs(kk,4:6)=[xy_curr_ids_locs(kk,1)+dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
        %Jy currents
        ind=num_nonair_cube+kk;
        all_panel_locs(ind,1:3)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)-dx*0.5 xy_curr_ids_locs(ind,3)];
        all_panel_locs(ind,4:6)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)+dx*0.5 xy_curr_ids_locs(ind,3)];
    end
    if(fl_profile == 1); disp(['Time for finding locations of surface panels ::: ', num2str(toc)]); end
    
    % finding unique panels
    tic
    tola=1e-12;
    unique_panel_coors_x=zeros(2*num_nonair_cube,3);
    unique_panel_coors_y=zeros(2*num_nonair_cube,3);
    dum_x=1;
    dum_y=1;
    for kk=1:num_nonair_cube
        ind=num_nonair_cube+kk;
        if (kk == 1)
            unique_panel_coors_x(dum_x,1:3)=all_panel_locs(kk,1:3);
            dum_x=dum_x+1;
            unique_panel_coors_x(dum_x,1:3)=all_panel_locs(kk,4:6);
            
            unique_panel_coors_y(dum_y,1:3)=all_panel_locs(ind,1:3);
            dum_y=dum_y+1;
            unique_panel_coors_y(dum_y,1:3)=all_panel_locs(ind,4:6);
            
        else
            %For enclosing Jx currents
            fl_found=0;
            %for ll=1:dum_x
            for ll=dum_x:-1:1
                if (abs(unique_panel_coors_x(ll,1)-all_panel_locs(kk,1))<tola && ...
                        abs(unique_panel_coors_x(ll,2)-all_panel_locs(kk,2))<tola && ...
                        abs(unique_panel_coors_x(ll,3)-all_panel_locs(kk,3))<tola)
                    fl_found=1;
                    break
                end
            end
            if (fl_found == 0)
                dum_x=dum_x+1;
                unique_panel_coors_x(dum_x,1:3)=all_panel_locs(kk,1:3);
            end
            
            fl_found=0;
            %for ll=1:dum_x
            for ll=dum_x:-1:1
                if (abs(unique_panel_coors_x(ll,1)-all_panel_locs(kk,4))<tola && ...
                        abs(unique_panel_coors_x(ll,2)-all_panel_locs(kk,5))<tola && ...
                        abs(unique_panel_coors_x(ll,3)-all_panel_locs(kk,6))<tola)
                    fl_found=1;
                    break
                end
            end
            if (fl_found == 0)
                dum_x=dum_x+1;
                unique_panel_coors_x(dum_x,1:3)=all_panel_locs(kk,4:6);
            end
            
            %For enclosing Jy currents
            
            fl_found=0;
            %for ll=1:dum_y
            for ll=dum_y:-1:1
                if (abs(unique_panel_coors_y(ll,1)-all_panel_locs(ind,1))<tola && ...
                        abs(unique_panel_coors_y(ll,2)-all_panel_locs(ind,2))<tola && ...
                        abs(unique_panel_coors_y(ll,3)-all_panel_locs(ind,3))<tola)
                    fl_found=1;
                    break
                end
            end
            if (fl_found == 0)
                dum_y=dum_y+1;
                unique_panel_coors_y(dum_y,1:3)=all_panel_locs(ind,1:3);
            end
            
            fl_found=0;
            %for ll=1:dum_y
            for ll=dum_y:-1:1
                if (abs(unique_panel_coors_y(ll,1)-all_panel_locs(ind,4))<tola && ...
                        abs(unique_panel_coors_y(ll,2)-all_panel_locs(ind,5))<tola && ...
                        abs(unique_panel_coors_y(ll,3)-all_panel_locs(ind,6))<tola)
                    fl_found=1;
                    break
                end
            end
            if (fl_found == 0)
                dum_y=dum_y+1;
                unique_panel_coors_y(dum_y,1:3)=all_panel_locs(ind,4:6);
            end
            
        end
    end
    if(fl_profile == 1); disp(['Time for finding unique surface panels ::: ', num2str(toc)]); end
    unique_panel_coors_x_dum=unique_panel_coors_x(1:dum_x,1:3);
    unique_panel_coors_y_dum=unique_panel_coors_y(1:dum_y,1:3);
    clear unique_panel_coors_x unique_panel_coors_y
    unique_panel_coors_x=unique_panel_coors_x_dum;
    unique_panel_coors_y=unique_panel_coors_y_dum;
    clear unique_panel_coors_x_dum unique_panel_coors_y_dum
    
    % Sorting unique panel coordinates in ascending order. Sort is performed
    % due to x coordinate (col=1), then if the values are the same at the first
    % column, then it sorts due to col=2 and so on...
    % The command sortrows(unique_panel_coors_x) is not working so brute force sorting!!!
    tic
    % sorting unique panels enclosing x-directed currents
    unique_panel_coors_x_tmp=unique_panel_coors_x;
    local_min_x=min(unique_panel_coors_x_tmp);
    local_panel_ids_x=zeros(size(unique_panel_coors_x_tmp,1),3);
    
    for kk=1:size(unique_panel_coors_x_tmp,1)
        local_panel_ids_x(kk,1) = (unique_panel_coors_x_tmp(kk,1) - local_min_x(1))/dx;
        local_panel_ids_x(kk,2) = (unique_panel_coors_x_tmp(kk,2) - local_min_x(2))/dx;
        local_panel_ids_x(kk,3) = (unique_panel_coors_x_tmp(kk,3) - local_min_x(3))/dx;
    end
    local_panel_ids_x = local_panel_ids_x+1; %convert to index format
    local_max_x=round(max(local_panel_ids_x));
    
    tmp_tens=zeros(local_max_x(1),local_max_x(2),local_max_x(3),3);
    tmp_tens_fl=zeros(local_max_x(1),local_max_x(2),local_max_x(3));
    % for finding the number of grids along all directions
    tola=1e-12;
    for kk=1:local_max_x(1) %x
        for ll=1:local_max_x(2) %y
            for mm=1:local_max_x(3) %z
                for nn=1:size(local_panel_ids_x,1)
                    if (abs(local_panel_ids_x(nn,1)-kk)<tola && ...
                            abs(local_panel_ids_x(nn,2)-ll)<tola && ...
                            abs(local_panel_ids_x(nn,3)-mm)<tola)
                        tmp_tens(kk,ll,mm,1:3)=unique_panel_coors_x_tmp(nn,1:3);
                        tmp_tens_fl(kk,ll,mm)=1;
                        break
                    end
                end
            end
        end
    end
    
    unique_panel_coors_x_tmp=0.0;
    dum=1;tola=1e-12;
    % Deposit tensor
    for kk=1:local_max_x(1) %x
        for ll=1:local_max_x(2) %y
            for mm=1:local_max_x(3) %z
                if (tmp_tens_fl(kk,ll,mm)==1)
                    unique_panel_coors_x_tmp(dum,1:3)=tmp_tens(kk,ll,mm,1:3);
                    dum=dum+1;
                end
            end
        end
    end
    
    if (size(unique_panel_coors_x_tmp,1) ~= size(unique_panel_coors_x,1))
        error('Something is fishy! the size of sorted array is not equal to non-sorted one')
    end
    unique_panel_coors_x=unique_panel_coors_x_tmp;
    
    % sorting unique panels enclosing y-directed currents
    unique_panel_coors_y_tmp=unique_panel_coors_y;
    local_min_y=min(unique_panel_coors_y_tmp);
    local_panel_ids_y=zeros(size(unique_panel_coors_y_tmp,1),3);
    
    
    for kk=1:size(unique_panel_coors_y_tmp,1)
        local_panel_ids_y(kk,1) = (unique_panel_coors_y_tmp(kk,1) - local_min_y(1))/dx;
        local_panel_ids_y(kk,2) = (unique_panel_coors_y_tmp(kk,2) - local_min_y(2))/dx;
        local_panel_ids_y(kk,3) = (unique_panel_coors_y_tmp(kk,3) - local_min_y(3))/dx;
    end
    local_panel_ids_y = local_panel_ids_y+1; %convert to index format
    local_max_y=round(max(local_panel_ids_y));
    
    tmp_tens=zeros(local_max_y(1),local_max_y(2),local_max_y(3),3);
    tmp_tens_fl=zeros(local_max_y(1),local_max_y(2),local_max_y(3));
    % for finding the number of grids along all directions
    tola=1e-12;
    for kk=1:local_max_y(1) %x
        for ll=1:local_max_y(2) %y
            for mm=1:local_max_y(3) %z
                for nn=1:size(local_panel_ids_y,1)
                    if (abs(local_panel_ids_y(nn,1)-kk)<tola && ...
                            abs(local_panel_ids_y(nn,2)-ll)<tola && ...
                            abs(local_panel_ids_y(nn,3)-mm)<tola)
                        tmp_tens(kk,ll,mm,1:3)=unique_panel_coors_y_tmp(nn,1:3);
                        tmp_tens_fl(kk,ll,mm)=1;
                        break
                    end
                end
            end
        end
    end
    
    unique_panel_coors_y_tmp=0.0;
    dum=1;tola=1e-12;
    % Deposit tensor
    for kk=1:local_max_y(1) %x
        for ll=1:local_max_y(2) %y
            for mm=1:local_max_y(3) %z
                if (tmp_tens_fl(kk,ll,mm)==1)
                    unique_panel_coors_y_tmp(dum,1:3)=tmp_tens(kk,ll,mm,1:3);
                    dum=dum+1;
                end
            end
        end
    end
    if (size(unique_panel_coors_y_tmp,1) ~= size(unique_panel_coors_y,1))
        error('Something is fishy! the size of sorted array is not equal to non-sorted one')
    end
    unique_panel_coors_y=unique_panel_coors_y_tmp;
    
    if(fl_profile == 1); disp(['Time for sorting panels ::: ', num2str(toc)]); end
    
    % Finding the ids of unique panels
    tic
    all_panels_ids=zeros(2*num_nonair_cube,2);
    for kk=1:num_nonair_cube
        
        for ll=1:size(unique_panel_coors_x,1)
            if (abs(unique_panel_coors_x(ll,1)-all_panel_locs(kk,1)) < tola && ...
                    abs(unique_panel_coors_x(ll,2)-all_panel_locs(kk,2)) < tola && ...
                    abs(unique_panel_coors_x(ll,3)-all_panel_locs(kk,3)) < tola)
                all_panels_ids(kk,1) = ll;
                break
            end
        end
        
        for ll=1:size(unique_panel_coors_x,1)
            if (abs(unique_panel_coors_x(ll,1)-all_panel_locs(kk,4)) < tola && ...
                    abs(unique_panel_coors_x(ll,2)-all_panel_locs(kk,5)) < tola && ...
                    abs(unique_panel_coors_x(ll,3)-all_panel_locs(kk,6)) < tola)
                all_panels_ids(kk,2) = ll;
                break
            end
        end
        
        ind=num_nonair_cube+kk;
        
        for ll=1:size(unique_panel_coors_y,1)
            if (abs(unique_panel_coors_y(ll,1)-all_panel_locs(ind,1)) < tola && ...
                    abs(unique_panel_coors_y(ll,2)-all_panel_locs(ind,2)) < tola && ...
                    abs(unique_panel_coors_y(ll,3)-all_panel_locs(ind,3)) < tola)
                all_panels_ids(ind,1) = size(unique_panel_coors_x,1)+ll;
                break
            end
        end
        
        for ll=1:size(unique_panel_coors_y,1)
            if (abs(unique_panel_coors_y(ll,1)-all_panel_locs(ind,4)) < tola && ...
                    abs(unique_panel_coors_y(ll,2)-all_panel_locs(ind,5)) < tola && ...
                    abs(unique_panel_coors_y(ll,3)-all_panel_locs(ind,6)) < tola)
                all_panels_ids(ind,2) = size(unique_panel_coors_x,1)+ll;
                break
            end
        end
    end
    
    if(fl_profile == 1); disp(['Time for finding the ids of unique panels ::: ', num2str(toc)]); end
    
    tic
    
    num_node_enc_x=size(unique_panel_coors_x,1); % # of panels enclosing Jx
    num_node_enc_y=size(unique_panel_coors_y,1); % # of panels enclosing Jy
    num_nodes=num_node_enc_x+num_node_enc_y;
    
    Ae = spalloc(num_nodes,3*num_nonair_cube,8*num_nonair_cube);
    
    infomem1 = whos('Ae');
    memestimated = (infomem1.bytes)/(1024*1024);
    disp(['Memory for Ae matrix (MB)::' , num2str(memestimated)]);
    
    % enter +1, leaves -1; first entry for leaving, second entry for entering
    for kk=1:2*num_nonair_cube
        Ae(all_panels_ids(kk,1),kk)=-1;
        Ae(all_panels_ids(kk,2),kk)=1;
    end
    %const_lin=dx/2;
    const_lin=1/2;
    Ae(1:num_node_enc_x,2*num_nonair_cube+1:3*num_nonair_cube) = ...
        abs(Ae(1:num_node_enc_x,1:num_nonair_cube)) * (const_lin); % for x comp of Jd
    
    Ae(num_node_enc_x+1:num_node_enc_x+num_node_enc_y,2*num_nonair_cube+1:3*num_nonair_cube) = ...
        abs(Ae(num_node_enc_x+1:num_node_enc_x+num_node_enc_y,num_nonair_cube+1:2*num_nonair_cube)) * (-const_lin); % for y comp of Jd
    
    if(fl_profile == 1); disp(['Time for generating sparse Ae matrix ::: ', num2str(toc)]); end;
    
else
    error('Invalid selection for Ae matrix generation: Set it to 0, 1, or 2...')
end

if(fl_profile == 1); disp(['Total time for generating Ae matrix ::: ', num2str(toc(tstart))]); end;

%% Visualize panel IDS with the number of unknown current it is enclosing
tic

if (fl_vis_only_voxels_with_ids == 1)
    figure; set(gca,'FontSize',24);
    xd = grid_intcon(:,:,:,1);yd = grid_intcon(:,:,:,2);zd = grid_intcon(:,:,:,3);
    h=plot3(xd(:), yd(:), zd(:), 'r*');
    set(h,'MarkerSize',10); xlabel('x'); ylabel('y'); zlabel('z');set(gca,'FontSize',24);
    hold on
    dum=1;
    for mm=1:N
        for ll=1:M
            for kk=1:L
                if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                        grid_intcon(kk,ll,mm,3) ~=0 )
                    coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                    if (mod(dum,freq_vis) == 0 || dum == 1 || dum == num_nonair_cube || kk == 1)
                        text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(dum));
                        %text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(num_nonair_cube+dum));
                    end
                    dum=dum+1;
                end
            end
        end
    end
    grid on
    axis tight
end


if (fl_vis_panel_locs_ids_x == 1)
    % panels enclosing Jx currents
    xd = grid_intcon(:,:,:,1);yd = grid_intcon(:,:,:,2);zd = grid_intcon(:,:,:,3);
    figure;
    set(gca,'FontSize',24);
    h=plot3(xd(:), yd(:), zd(:), 'r*');
    set(h,'MarkerSize',10); xlabel('x'); ylabel('y'); zlabel('z');
    
    %     if (isempty(region_focus) == 0) % no focus
    %         % plot the intervals of focus region
    %         hold on
    %         h=plot3(region_focus(1,1), region_focus(2,1), region_focus(3,1), 'ko');
    %         set(h,'MarkerSize',10);
    %         hold on
    %         h=plot3(region_focus(1,2), region_focus(2,2), region_focus(3,2), 'ko');
    %         set(h,'MarkerSize',10);
    %         hold on
    %         h=plot3(exc_pnt(1), exc_pnt(2), exc_pnt(3), 'k>');
    %         set(h,'MarkerSize',10);
    %     end
    tola=1e-12;
    dum=1;
    for mm=1:N % Attention:: Ordering and numbering were changed!!!
        for ll=1:M
            for kk=1:L
                if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                        grid_intcon(kk,ll,mm,3) ~=0 )
                    hold on
                    h=plot3(all_panel_locs(dum,1),all_panel_locs(dum,2),all_panel_locs(dum,3),'b.');
                    set(h,'MarkerSize',10);
                    hold on
                    h=plot3(all_panel_locs(dum,4),all_panel_locs(dum,5),all_panel_locs(dum,6),'b.');
                    set(h,'MarkerSize',10);
                    hold on
                    %                     if (isempty(region_focus) == 1) % no focus
                    if (mod(dum,freq_vis) == 0 || dum == 1 || dum == num_nonair_cube || kk == 1)
                        coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                        h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(dum));set(h,'FontSize',24);
                        h=text(all_panel_locs(dum,1),all_panel_locs(dum,2),all_panel_locs(dum,3),num2str(all_panels_ids(dum,1)));set(h,'FontSize',24);
                        h=text(all_panel_locs(dum,4),all_panel_locs(dum,5),all_panel_locs(dum,6),num2str(all_panels_ids(dum,2)));set(h,'FontSize',24);
                    end
                    %                     else
                    %                         if (all_panel_locs(dum,1) >= region_focus(1,1)-tola && all_panel_locs(dum,1) <= region_focus(1,2)+tola && ...
                    %                                 all_panel_locs(dum,2) >= region_focus(2,1)-tola && all_panel_locs(dum,2) <= region_focus(2,2)+tola && ...
                    %                                 all_panel_locs(dum,3) >= region_focus(3,1)-tola && all_panel_locs(dum,3) <= region_focus(3,2)+tola)
                    %                             coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                    %                             h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(dum));set(h,'FontSize',24);
                    %                             h=text(all_panel_locs(dum,1),all_panel_locs(dum,2),all_panel_locs(dum,3),num2str(all_panels_ids(dum,1)));set(h,'FontSize',24);
                    %                             h=text(all_panel_locs(dum,4),all_panel_locs(dum,5),all_panel_locs(dum,6),num2str(all_panels_ids(dum,2)));set(h,'FontSize',24);
                    %                             if (isempty(exc_pnt) == 0)
                    %                                 coors_search=[coors_search;[all_panels_ids(dum,1) all_panel_locs(dum,1) all_panel_locs(dum,2) all_panel_locs(dum,3)]];
                    %                                 coors_search=[coors_search;[all_panels_ids(dum,2) all_panel_locs(dum,4) all_panel_locs(dum,5) all_panel_locs(dum,6)]];
                    %                             end
                    %                         end
                    %                     end
                    
                    dum=dum+1;
                end
            end
        end
    end
    
    title('IDs of nodes enclosing Jx currents');
    grid on; set(gca,'FontSize',24);
    %view(0,0)
    view(-25,45)
    %     if (isempty(region_focus) == 0)
    %         xlim([region_focus(1,1) region_focus(1,2)])
    %         ylim([region_focus(2,1) region_focus(2,2)])
    %         zlim([region_focus(3,1) region_focus(3,2)])
    %     end
    %     nodes_exc=zeros(1,2);
    %     if (isempty(exc_pnt) == 0)
    %         dists=zeros(size(coors_search,1),1);
    %         for kk=1:size(coors_search,1)
    %             dists(kk)=norm(exc_pnt-coors_search(kk,2:4));
    %         end
    %         [Y,ind_vect] = sort(dists);
    %         %ind_vect(1:4)
    %         %coors_search
    %         tmp_inds=coors_search(ind_vect(1:4),1)
    %         if( fl_meth_comp==0 || fl_meth_comp==1 )
    %             nodes_exc(1)=min(tmp_inds);
    %             nodes_exc(2)=max(tmp_inds)-1;
    %         elseif (fl_meth_comp==2 )
    %             %nodes_exc(1)=max(tmp_inds);
    %             %nodes_exc(2)=max(tmp_inds)-1;
    %             nodes_exc(1)=tmp_inds(1);
    %             nodes_exc(2)=tmp_inds(2);
    %
    %         end
    %         disp(['Ids of nodes for connecting current source ::: ',num2str(nodes_exc)])
    %     end
end

if (fl_vis_panel_locs_ids_y == 1)
    % panels enclosing Jy currents
    xd = grid_intcon(:,:,:,1);yd = grid_intcon(:,:,:,2);zd = grid_intcon(:,:,:,3);
    figure;
    set(gca,'FontSize',24);
    h=plot3(xd(:), yd(:), zd(:), 'r*');
    set(h,'MarkerSize',10); xlabel('x'); ylabel('y'); zlabel('z');
    dum=1;
    for mm=1:N % Attention:: Ordering and numbering were changed!!!
        for ll=1:M
            for kk=1:L
                if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                        grid_intcon(kk,ll,mm,3) ~=0 )
                    hold on
                    ind=num_nonair_cube+dum;
                    
                    h=plot3(all_panel_locs(ind,1),all_panel_locs(ind,2),all_panel_locs(ind,3),'b.');
                    set(h,'MarkerSize',10);
                    hold on
                    h=plot3(all_panel_locs(ind,4),all_panel_locs(ind,5),all_panel_locs(ind,6),'b.');
                    set(h,'MarkerSize',10);
                    hold on
                    
                    if (mod(dum,freq_vis) == 0 || dum == 1 || dum == num_nonair_cube || kk == 1)
                        coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                        h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(ind));set(h,'FontSize',24);
                        h=text(all_panel_locs(ind,1),all_panel_locs(ind,2),all_panel_locs(ind,3),num2str(all_panels_ids(ind,1)));set(h,'FontSize',24);
                        h=text(all_panel_locs(ind,4),all_panel_locs(ind,5),all_panel_locs(ind,6),num2str(all_panels_ids(ind,2)));set(h,'FontSize',24);
                    end
                    
                    dum=dum+1;
                end
            end
        end
    end
    title('IDs of nodes enclosing Jy currents');
    grid on; set(gca,'FontSize',24);
    view(-25,45)
end


if (fl_vis_panel_locs_ids_z == 1)
    % panels enclosing Jz currents
    xd = grid_intcon(:,:,:,1);yd = grid_intcon(:,:,:,2);zd = grid_intcon(:,:,:,3);
    figure;
    set(gca,'FontSize',24);
    h=plot3(xd(:), yd(:), zd(:), 'r*');
    set(h,'MarkerSize',10); xlabel('x'); ylabel('y'); zlabel('z');
    dum=1;
    for mm=1:N % Attention:: Ordering and numbering were changed!!!
        for ll=1:M
            for kk=1:L
                if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                        grid_intcon(kk,ll,mm,3) ~=0 )
                    hold on
                    ind=2*num_nonair_cube+dum;
                    
                    h=plot3(all_panel_locs(ind,1),all_panel_locs(ind,2),all_panel_locs(ind,3),'b.');
                    set(h,'MarkerSize',10);
                    hold on
                    h=plot3(all_panel_locs(ind,4),all_panel_locs(ind,5),all_panel_locs(ind,6),'b.');
                    set(h,'MarkerSize',10);
                    hold on
                    
                    if (mod(dum,freq_vis) == 0 || dum == 1 || dum == num_nonair_cube || kk == 1)
                        coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                        h=text(coor_tmp(1),coor_tmp(2),coor_tmp(3),num2str(ind));set(h,'FontSize',24);
                        h=text(all_panel_locs(ind,1),all_panel_locs(ind,2),all_panel_locs(ind,3),num2str(all_panels_ids(ind,1)));set(h,'FontSize',24);
                        h=text(all_panel_locs(ind,4),all_panel_locs(ind,5),all_panel_locs(ind,6),num2str(all_panels_ids(ind,2)));set(h,'FontSize',24);
                    end
                    
                    dum=dum+1;
                end
            end
        end
    end
    title('IDs of nodes enclosing Jz currents');
    grid on; set(gca,'FontSize',24);
    %view(0,0)
    view(-25,45)
end
if(fl_profile == 1); disp(['Time for visualizing panels and currents ::: ', num2str(toc)]); end

%% Finding the nodes on which excitation and ground points are defined

disp('-----------------------------------------------------')
disp('Finding the IDs of port nodes... ')


if(fl_vis_exc_grnd == 1)
    if (isempty(pnt_exc) == 0 || isempty(pnt_grnd) == 0)
        figure;
    end
end

if (isempty(pnt_exc) == 0)
    tic
    
    % finding the ids of non-empty voxels
    
    voxel_ids=zeros(N,M,L);
    dum=1;
    for mm=1:N % Attention:: Ordering and numbering were changed!!!
        for ll=1:M
            for kk=1:L
                if ( grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                        grid_intcon(kk,ll,mm,3) ~=0 )
                   
                    voxel_ids(mm,ll,kk)=dum;
                    dum=dum+1;
                    
                end
            end
        end
    end
    
    
    
    tola=1e-12;
    nodeid_4_injectcurr=zeros(size(pnt_exc,1),1);
    %unkid_4_injectcurr=zeros(size(pnt_exc,1),1);
    for nn=1:size(pnt_exc,1)
        pnt_coor = pnt_exc(nn,:);
        dum=1;
        pnt_found=0;
        
        tent_x=ceil((pnt_coor(1)-grid_intcon(1,1,1,1))/dx);
        tent_y=ceil((pnt_coor(2)-grid_intcon(1,1,1,2))/dx);
        tent_z=ceil((pnt_coor(3)-grid_intcon(1,1,1,3))/dx);
        
        x_beg=tent_x-5; x_end=tent_x+5;
        y_beg=tent_y-5; y_end=tent_y+5;
        z_beg=tent_z-5; z_end=tent_z+5;
        
        if (x_beg>L); x_beg=L; end; if (x_beg<1); x_beg=1; end;
        if (x_end>L); x_end=L; end; if (x_end<1); x_end=1; end;
        
        if (y_beg>M); y_beg=M; end; if (y_beg<1); y_beg=1; end;
        if (y_end>M); y_end=M; end; if (y_end<1); y_end=1; end;
        
        if (z_beg>N); z_beg=N; end; if (z_beg<1); z_beg=1; end;
        if (z_end>N); z_end=N; end; if (z_end<1); z_end=1; end;

        for mm=z_beg:z_end % Attention:: Ordering and numbering were changed!!!
            for ll=y_beg:y_end
                for kk=x_beg:x_end
        %for mm=1:N % Attention:: Ordering and numbering were changed!!!
        %    for ll=1:M
        %        for kk=1:L
                    if ( grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                            grid_intcon(kk,ll,mm,3) ~=0 )
                        
                        % starts here ~~~~~
                        dum=squeeze(voxel_ids(mm,ll,kk));
                        
                        [pnt_found,pnt_id] = lse_find_nodeid(all_panel_locs,all_panels_ids,pnt_coor,num_nonair_cube,tola,fl_vis_exc_grnd,dum);
                        
                        % ends here ~~~~~
                        
                        if(pnt_found == 1)
                            nodeid_4_injectcurr(nn) = pnt_id;
                            break
                        end
                        %dum=dum+1;
                        
                    end
                end
                if(pnt_found == 1)
                    break
                end
            end
            if(pnt_found == 1)
                break
            end
        end
        if (nodeid_4_injectcurr(nn) == 0)
            error(['node id for exc pnt ::: ',num2str(pnt_coor),' could not find!!!'])
        end
        
    end
    if(fl_vis_exc_grnd == 1)
        xlabel('x');ylabel('y');zlabel('z');set(gca,'FontSize',24);
        axis tight; grid on; view(40,20)
    end
    if(fl_profile == 1); disp(['Time for finding & visualizing excitation nodes ::: ', num2str(toc)]); end
end



if (isempty(pnt_grnd) == 0)
    tic
    tola=1e-12;
    nodeid_4_grnd=zeros(size(pnt_grnd,1),1);
    %unkid_4_grnd=zeros(size(pnt_grnd,1),1);
    for nn=1:size(pnt_grnd,1)
        pnt_coor = pnt_grnd(nn,:);
        dum=1;
        pnt_found=0;
        
        tent_x=ceil((pnt_coor(1)-grid_intcon(1,1,1,1))/dx);
        tent_y=ceil((pnt_coor(2)-grid_intcon(1,1,1,2))/dx);
        tent_z=ceil((pnt_coor(3)-grid_intcon(1,1,1,3))/dx);
        
        x_beg=tent_x-5; x_end=tent_x+5;
        y_beg=tent_y-5; y_end=tent_y+5;
        z_beg=tent_z-5; z_end=tent_z+5;
        
        if (x_beg>L); x_beg=L; end; if (x_beg<1); x_beg=1; end;
        if (x_end>L); x_end=L; end; if (x_end<1); x_end=1; end;
        
        if (y_beg>M); y_beg=M; end; if (y_beg<1); y_beg=1; end;
        if (y_end>M); y_end=M; end; if (y_end<1); y_end=1; end;
        
        if (z_beg>N); z_beg=N; end; if (z_beg<1); z_beg=1; end;
        if (z_end>N); z_end=N; end; if (z_end<1); z_end=1; end;
        
        for mm=z_beg:z_end % Attention:: Ordering and numbering were changed!!!
            for ll=y_beg:y_end
                for kk=x_beg:x_end
        
%         for mm=1:N % Attention:: Ordering and numbering were changed!!!
%             for ll=1:M
%                 for kk=1:L
                    if ( grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                            grid_intcon(kk,ll,mm,3) ~=0 )
                        
                        % starts here ~~~~~
                        dum=squeeze(voxel_ids(mm,ll,kk));
                        [pnt_found,pnt_id] = lse_find_nodeid(all_panel_locs,all_panels_ids,pnt_coor,num_nonair_cube,tola,fl_vis_exc_grnd,dum);
                        
                        % ends here ~~~~~
                        
                        if(pnt_found == 1)
                            nodeid_4_grnd(nn) = pnt_id;
                            break
                        end
                        
                        %dum=dum+1;
                        
                    end
                end
                if(pnt_found == 1)
                    break
                end
            end
            if(pnt_found == 1)
                break
            end
        end
        if (nodeid_4_grnd(nn) == 0)
            error(['node id for ground pnt ::: ',num2str(pnt_coor),' could not find!!!'])
        end
        if(fl_vis_exc_grnd == 1)
            xlabel('x');ylabel('y');zlabel('z');set(gca,'FontSize',24);
            axis tight; grid on; view(40,20)
        end
    end
    if(fl_profile == 1); disp(['Time for finding & visualizing ground nodes ::: ', num2str(toc)]); end
end

if (isempty(pnt_well_cond) == 0)
    tic
    tola=1e-12;
    nodeid_4_well_cond=zeros(size(pnt_well_cond,1),1);
    for nn=1:size(pnt_well_cond,1)
        pnt_coor = pnt_well_cond(nn,:);
        dum=1;
        pnt_found=0;
        
        tent_x=ceil((pnt_coor(1)-grid_intcon(1,1,1,1))/dx);
        tent_y=ceil((pnt_coor(2)-grid_intcon(1,1,1,2))/dx);
        tent_z=ceil((pnt_coor(3)-grid_intcon(1,1,1,3))/dx);
        
        x_beg=tent_x-5; x_end=tent_x+5;
        y_beg=tent_y-5; y_end=tent_y+5;
        z_beg=tent_z-5; z_end=tent_z+5;
        
        if (x_beg>L); x_beg=L; end; if (x_beg<1); x_beg=1; end;
        if (x_end>L); x_end=L; end; if (x_end<1); x_end=1; end;
        
        if (y_beg>M); y_beg=M; end; if (y_beg<1); y_beg=1; end;
        if (y_end>M); y_end=M; end; if (y_end<1); y_end=1; end;
        
        if (z_beg>N); z_beg=N; end; if (z_beg<1); z_beg=1; end;
        if (z_end>N); z_end=N; end; if (z_end<1); z_end=1; end;
        
        for mm=z_beg:z_end % Attention:: Ordering and numbering were changed!!!
            for ll=y_beg:y_end
                for kk=x_beg:x_end
        
        %for mm=1:N % Attention:: Ordering and numbering were changed!!!
        %    for ll=1:M
        %        for kk=1:L
                    if ( grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                            grid_intcon(kk,ll,mm,3) ~=0 )
                        
                        % starts here ~~~~~
                        dum=squeeze(voxel_ids(mm,ll,kk));
                        [pnt_found,pnt_id] = lse_find_nodeid(all_panel_locs,all_panels_ids,pnt_coor,num_nonair_cube,tola,fl_vis_exc_grnd,dum);
                        
                        % ends here ~~~~~
                        
                        if(pnt_found == 1)
                            nodeid_4_well_cond(nn) = pnt_id;
                            break
                        end
                        
                        %dum=dum+1;
                    end
                end
                if(pnt_found == 1)
                    break
                end
            end
            if(pnt_found == 1)
                break
            end
        end
        if (nodeid_4_well_cond(nn) == 0)
            error(['node id for well-conditioner ground pnt ::: ',num2str(pnt_coor),' could not find!!!'])
        end
        if(fl_vis_exc_grnd == 1)
            xlabel('x');ylabel('y');zlabel('z');set(gca,'FontSize',24);
            axis tight; grid on; view(40,20)
        end
    end
    if(fl_profile == 1); disp(['Time for finding & visualizing ground nodes ::: ', num2str(toc)]); end
end

% putting the node ids of ports nodes in structure

% puting node/current ids of each port in a structure
nodeid_lft=cell(num_ports,1);
nodeid_rght=cell(num_ports,1);

dum_exc=0;
dum_grnd=0;
for kk=1:num_ports
    nodeid_lft{kk}(:) = nodeid_4_injectcurr(dum_exc+1:dum_exc+num_node_each_port(kk,2));
    dum_exc=dum_exc+num_node_each_port(kk,2);
    
    nodeid_rght{kk}(:) = nodeid_4_grnd(dum_grnd+1:dum_grnd+num_node_each_port(kk,1));
    dum_grnd=dum_grnd+num_node_each_port(kk,1);
end

nodeid_wlcond = nodeid_4_well_cond;

disp('Done... Finding the IDs of port nodes... ')
disp('-----------------------------------------------------')

%if (fl_check_ports == 1); error('Just plotting port nodes... No Simulation...'); end

num_node = size(Ae,1);
num_curr = size(Ae,2);
disp('-----------------------------------------------------')
disp(['Number of current unknowns ::: ',num2str(num_curr)])
disp(['Number of potential unknowns ::: ',num2str(num_node)])
disp(['Number of total unknowns ::: ',num2str(num_curr+num_node)])
disp('-----------------------------------------------------')

