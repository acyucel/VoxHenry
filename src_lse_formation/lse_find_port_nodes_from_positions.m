function [nodeid_lft, nodeid_rght, nodeid_wlcond] = lse_find_port_nodes_from_positions(all_panel_locs, all_panels_ids, grid_intcon, L, M, N, dx, pnt_lft, pnt_rght, pnt_well_cond, fl_check_ports)
% Finding the nodes on which excitation and ground points are defined,
% based on the physical location in space of the port nodes

fl_profile = 0;

% Inputs for visualization
fl_vis_exc_grnd = fl_check_ports;

nodeid_4_injectcurr=[];
nodeid_4_grnd=[];
nodeid_4_well_cond=[];


num_nonair_cube = size(all_panels_ids,1);

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
    nodeid_4_injectcurr=zeros(size(pnt_exc,1),3);
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
                        
                        [pnt_found,pnt_id,side] = lse_find_nodeid(all_panel_locs,all_panels_ids,pnt_coor,num_nonair_cube,tola,fl_vis_exc_grnd,dum);
                        
                        % ends here ~~~~~
                        
                        if(pnt_found == 1)
                            nodeid_4_injectcurr(nn,:) = [pnt_id dum side];
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
        if (nodeid_4_injectcurr(nn,1) == 0)
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
    nodeid_4_grnd=zeros(size(pnt_grnd,1),3);
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
                        [pnt_found,pnt_id,side] = lse_find_nodeid(all_panel_locs,all_panels_ids,pnt_coor,num_nonair_cube,tola,fl_vis_exc_grnd,dum);
                        
                        % ends here ~~~~~
                        
                        if(pnt_found == 1)
                            nodeid_4_grnd(nn,:) = [pnt_id dum side];
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
        if (nodeid_4_grnd(nn,1) == 0)
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
    nodeid_4_well_cond=zeros(size(pnt_well_cond,1),3);
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
                        [pnt_found,pnt_id,side] = lse_find_nodeid(all_panel_locs,all_panels_ids,pnt_coor,num_nonair_cube,tola,fl_vis_exc_grnd,dum);
                        
                        % ends here ~~~~~
                        
                        if(pnt_found == 1)
                            nodeid_4_well_cond(nn,:) = [pnt_id dum side];
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
        if (nodeid_4_well_cond(nn,1) == 0)
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
    % changed from "nodeid_lft{kk}(:)" to "nodeid_lft{kk}" to work under Octave
    %nodeid_lft{kk}(:) = nodeid_4_injectcurr(dum_exc+1:dum_exc+num_node_each_port(kk,2));
    nodeid_lft{kk} = nodeid_4_injectcurr(dum_exc+1:dum_exc+num_node_each_port(kk,2),:);
    dum_exc=dum_exc+num_node_each_port(kk,2);
    
    %nodeid_rght{kk}(:) = nodeid_4_grnd(dum_grnd+1:dum_grnd+num_node_each_port(kk,1));
    nodeid_rght{kk} = nodeid_4_grnd(dum_grnd+1:dum_grnd+num_node_each_port(kk,1),:);
    dum_grnd=dum_grnd+num_node_each_port(kk,1);
end

nodeid_wlcond = nodeid_4_well_cond;

disp('Done... Finding the IDs of port nodes... ')
disp('-----------------------------------------------------')


