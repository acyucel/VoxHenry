function [all_panel_locs] = lse_find_node_locations(num_nonair_cube, grid_intcon, L, M, N, dx, all_panels_ids)
% find the position in space of the six nodes of each non-air voxel 

fl_profile = 0;


% Inputs for visualization
fl_vis_only_voxels_with_ids=0;
fl_vis_panel_locs_ids_x=0;
fl_vis_panel_locs_ids_y=0;
fl_vis_panel_locs_ids_z=0;

num_elem_allowed = 5000;
freq_vis = 1;
if (num_nonair_cube > num_elem_allowed)% visualization of boxes is slow when nD > num_elem_allowed
    freq_vis = 25; % every xx element unknown id will be written.
end


%
% finding locations of surface panels (coordinates)
% This is needed for debug visualization and for assigning the port nodes
% 

% centers of non-air voxels
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

% actual locations in space of the surface panels (actually, the node positions)
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




%
%% Visualize panel IDS with the number of unknown current it is enclosing
%

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
                        h=text(all_panel_locs(ind,1),all_panel_locs(ind,2),all_panel_locs(ind,3),num2str(all_panels_ids(dum,1)));set(h,'FontSize',24);
                        h=text(all_panel_locs(ind,4),all_panel_locs(ind,5),all_panel_locs(ind,6),num2str(all_panels_ids(dum,2)));set(h,'FontSize',24);
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
                        h=text(all_panel_locs(ind,1),all_panel_locs(ind,2),all_panel_locs(ind,3),num2str(all_panels_ids(dum,1)));set(h,'FontSize',24);
                        h=text(all_panel_locs(ind,4),all_panel_locs(ind,5),all_panel_locs(ind,6),num2str(all_panels_ids(dum,2)));set(h,'FontSize',24);
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






