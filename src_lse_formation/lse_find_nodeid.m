function [pnt_found,pnt_id,side] = lse_find_nodeid(all_panel_locs,all_panels_ids,pnt_coor,num_nonair_cube,tola,fl_vis_exc_grnd,dum)

pnt_found = 0;
pnt_id = 0;
side = 0;
if ( abs(all_panel_locs(dum,1)-pnt_coor(1)) < tola && ...
        abs(all_panel_locs(dum,2)-pnt_coor(2)) < tola && ...
        abs(all_panel_locs(dum,3)-pnt_coor(3)) < tola)
    
    side = 1;
    pnt_id=all_panels_ids(dum,side);

    % plot the point
    if(fl_vis_exc_grnd == 1)
        hold on
        h=plot3(all_panel_locs(dum,1),all_panel_locs(dum,2),all_panel_locs(dum,3),'r.');set(h,'MarkerSize',10);
        hold on
        h=text(all_panel_locs(dum,1),all_panel_locs(dum,2),all_panel_locs(dum,3),num2str(all_panels_ids(dum,1)));set(h,'FontSize',24);
    end
    pnt_found=1;
    return
    
elseif ( abs(all_panel_locs(dum,4)-pnt_coor(1)) < tola && ...
        abs(all_panel_locs(dum,5)-pnt_coor(2)) < tola && ...
        abs(all_panel_locs(dum,6)-pnt_coor(3)) < tola)
    
    side = 2;
    pnt_id=all_panels_ids(dum,side);

    % plot the point
    if(fl_vis_exc_grnd == 1)
        hold on
        h=plot3(all_panel_locs(dum,4),all_panel_locs(dum,5),all_panel_locs(dum,6),'r.');set(h,'MarkerSize',10);
        hold on
        h=text(all_panel_locs(dum,4),all_panel_locs(dum,5),all_panel_locs(dum,6),num2str(all_panels_ids(dum,2)));set(h,'FontSize',24);
    end
    
    pnt_found=1;
    if(pnt_found == 1)
        return
    end
    
    
elseif ( abs(all_panel_locs(num_nonair_cube+dum,1)-pnt_coor(1)) < tola && ...
        abs(all_panel_locs(num_nonair_cube+dum,2)-pnt_coor(2)) < tola && ...
        abs(all_panel_locs(num_nonair_cube+dum,3)-pnt_coor(3)) < tola)
    
    side = 3;
    pnt_id=all_panels_ids(dum,side);

    % plot the point
    if(fl_vis_exc_grnd == 1)
        hold on
        h=plot3(all_panel_locs(num_nonair_cube+dum,1),all_panel_locs(num_nonair_cube+dum,2),all_panel_locs(num_nonair_cube+dum,3),'r.');set(h,'MarkerSize',10);
        hold on
        h=text(all_panel_locs(num_nonair_cube+dum,1),all_panel_locs(num_nonair_cube+dum,2),all_panel_locs(num_nonair_cube+dum,3),num2str(all_panels_ids(dum,3)));set(h,'FontSize',24);
    end
    pnt_found=1;
    return
    
elseif ( abs(all_panel_locs(num_nonair_cube+dum,4)-pnt_coor(1)) < tola && ...
        abs(all_panel_locs(num_nonair_cube+dum,5)-pnt_coor(2)) < tola && ...
        abs(all_panel_locs(num_nonair_cube+dum,6)-pnt_coor(3)) < tola)
    
    side = 4;
    pnt_id=all_panels_ids(dum,side);

    % plot the point
    if(fl_vis_exc_grnd == 1)
        hold on
        h=plot3(all_panel_locs(num_nonair_cube+dum,4),all_panel_locs(num_nonair_cube+dum,5),all_panel_locs(num_nonair_cube+dum,6),'r.');set(h,'MarkerSize',10);
        hold on
        h=text(all_panel_locs(num_nonair_cube+dum,4),all_panel_locs(num_nonair_cube+dum,5),all_panel_locs(num_nonair_cube+dum,6),num2str(all_panels_ids(dum,4)));set(h,'FontSize',24);
    end
    
    pnt_found=1;
    if(pnt_found == 1)
        return
    end
    
elseif ( abs(all_panel_locs(2*num_nonair_cube+dum,1)-pnt_coor(1)) < tola && ...
        abs(all_panel_locs(2*num_nonair_cube+dum,2)-pnt_coor(2)) < tola && ...
        abs(all_panel_locs(2*num_nonair_cube+dum,3)-pnt_coor(3)) < tola)
    
    side = 5;
    pnt_id=all_panels_ids(dum,side);

    % plot the point
    if(fl_vis_exc_grnd == 1)
        hold on
        h=plot3(all_panel_locs(2*num_nonair_cube+dum,1),all_panel_locs(2*num_nonair_cube+dum,2),all_panel_locs(2*num_nonair_cube+dum,3),'r.');set(h,'MarkerSize',10);
        hold on
        h=text(all_panel_locs(2*num_nonair_cube+dum,1),all_panel_locs(2*num_nonair_cube+dum,2),all_panel_locs(2*num_nonair_cube+dum,3),num2str(all_panels_ids(dum,5)));set(h,'FontSize',24);
    end
    pnt_found=1;
    return
    
elseif ( abs(all_panel_locs(2*num_nonair_cube+dum,4)-pnt_coor(1)) < tola && ...
        abs(all_panel_locs(2*num_nonair_cube+dum,5)-pnt_coor(2)) < tola && ...
        abs(all_panel_locs(2*num_nonair_cube+dum,6)-pnt_coor(3)) < tola)
    
    side = 6;
    pnt_id=all_panels_ids(dum,side);

    % plot the point
    if(fl_vis_exc_grnd == 1)
        hold on
        h=plot3(all_panel_locs(2*num_nonair_cube+dum,4),all_panel_locs(2*num_nonair_cube+dum,5),all_panel_locs(2*num_nonair_cube+dum,6),'r.');set(h,'MarkerSize',10);
        hold on
        h=text(all_panel_locs(2*num_nonair_cube+dum,4),all_panel_locs(2*num_nonair_cube+dum,5),all_panel_locs(2*num_nonair_cube+dum,6),num2str(all_panels_ids(dum,6)));set(h,'FontSize',24);
    end
    
    pnt_found=1;
    if(pnt_found == 1)
        return
    end
    
end