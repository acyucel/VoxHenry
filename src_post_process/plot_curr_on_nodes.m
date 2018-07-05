function plot_curr_on_nodes(slct_plane,slct_cut,dx,nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned)

fl_ind_or_phys_len=0; % if 1 axes are physical length, otherwise index of node along that dimension
tola=1e-12;
fl_plot_surf_imagesc=2; % 1 for surf, 2 for imagesc

fl_slct_part=3; %1-real,2-imag,3-abs

if (fl_slct_part == 1) % real part of data for plotting
    nodes_w_currs_x_aligned(:,4)=real(nodes_w_currs_x_aligned(:,4));
    nodes_w_currs_y_aligned(:,4)=real(nodes_w_currs_y_aligned(:,4));
    nodes_w_currs_z_aligned(:,4)=real(nodes_w_currs_z_aligned(:,4));
    str_Jx='{\it{Re}}(J_x) ';str_Jy='{\it{Re}}(J_y) ';str_Jz='{\it{Re}}(J_z) ';
elseif (fl_slct_part == 2) % imag part of data for plotting
    nodes_w_currs_x_aligned(:,4)=imag(nodes_w_currs_x_aligned(:,4));
    nodes_w_currs_y_aligned(:,4)=imag(nodes_w_currs_y_aligned(:,4));
    nodes_w_currs_z_aligned(:,4)=imag(nodes_w_currs_z_aligned(:,4));
    str_Jx='{\it{Im}}(J_x) ';str_Jy='{\it{Im}}(J_y) ';str_Jz='{\it{Im}}(J_z) ';
elseif (fl_slct_part == 3)% abs value of data for plotting
    nodes_w_currs_x_aligned(:,4)=abs(nodes_w_currs_x_aligned(:,4));
    nodes_w_currs_y_aligned(:,4)=abs(nodes_w_currs_y_aligned(:,4));
    nodes_w_currs_z_aligned(:,4)=abs(nodes_w_currs_z_aligned(:,4));
    str_Jx='|J_x| ';str_Jy='|J_y| ';str_Jz='|J_z| ';
end

if (slct_plane == 'xy')
    eks_ch=3;
    eks_fix=[1 2];
elseif (slct_plane == 'xz')
    eks_ch=2;
    eks_fix=[1 3];
elseif (slct_plane == 'yz')
    eks_ch=1;
    eks_fix=[2 3];
end


figure

for comp=1:3 % for each component
    if (slct_plane == 'xy')
        % if we plot Jx and Jy, those are on the same plane defined by
        % select-cut
        % however if we plot Jz, that is on plane passing through slct_cut+dx/2
        if (comp == 1) % Jx
            slct_cut2=(slct_cut*dx)-dx/2;
            nodes_w_currs_dum=nodes_w_currs_x_aligned;
        elseif (comp ==2) % Jy
            slct_cut2=(slct_cut*dx)-dx/2;
            nodes_w_currs_dum=nodes_w_currs_y_aligned;
        elseif (comp == 3) %Jz
            slct_cut2=(slct_cut*dx);
            nodes_w_currs_dum=nodes_w_currs_z_aligned;
        end
    elseif (slct_plane == 'xz')
        if (comp == 1) % Jx
            slct_cut2=(slct_cut*dx)-dx/2;
            nodes_w_currs_dum=nodes_w_currs_x_aligned;
        elseif (comp ==2) % Jy
            slct_cut2=(slct_cut*dx);
            nodes_w_currs_dum=nodes_w_currs_y_aligned;
        elseif (comp == 3) %Jz
            slct_cut2=(slct_cut*dx)-dx/2;
            nodes_w_currs_dum=nodes_w_currs_z_aligned;
        end
        
    elseif (slct_plane == 'yz')
        if (comp == 1) % Jx
            slct_cut2=(slct_cut*dx);
            nodes_w_currs_dum=nodes_w_currs_x_aligned;
        elseif (comp ==2) % Jy
            slct_cut2=(slct_cut*dx)-dx/2;
            nodes_w_currs_dum=nodes_w_currs_y_aligned;
        elseif (comp == 3) %Jz
            slct_cut2=(slct_cut*dx)-dx/2;
            nodes_w_currs_dum=nodes_w_currs_z_aligned;
        end
    end
    

    
    tmp_book_keep=zeros(size(nodes_w_currs_dum,1),1);
    
    % find the z aligned nodes on specified cut
    for kk=1:size(nodes_w_currs_dum,1)
        if (abs(nodes_w_currs_dum(kk,eks_ch)-slct_cut2)<tola)
            tmp_book_keep(kk)=1;
        end
    end
    
    currs_plot= zeros(sum(tmp_book_keep),4);
    dum=1;
    for kk=1:size(nodes_w_currs_dum,1)
        if (tmp_book_keep(kk) == 1)
            currs_plot(dum,1:2)=nodes_w_currs_dum(kk,eks_fix);
            currs_plot(dum,3)=nodes_w_currs_dum(kk,4);
            dum=dum+1;
        end
    end
    
    % defining a grid on the cut and putting values in matrix
    num_smpl_x=round((max(currs_plot(:,1))-min(currs_plot(:,1)))/dx)+1;
    num_smpl_y=round((max(currs_plot(:,2))-min(currs_plot(:,2)))/dx)+1;
    
    xlin=linspace(min(currs_plot(:,1)),max(currs_plot(:,1)),num_smpl_x);
    ylin=linspace(min(currs_plot(:,2)),max(currs_plot(:,2)),num_smpl_y);
    
    [XX,YY]=meshgrid(xlin,ylin);
    org_x=min(currs_plot(:,1));
    org_y=min(currs_plot(:,2));
    curr_on_grid=zeros(size(XX,1),size(XX,2));
    for kk=1:size(currs_plot,1)
        indx=round((currs_plot(kk,1)-org_x)/dx)+1;
        indy=round((currs_plot(kk,2)-org_y)/dx)+1;
        curr_on_grid(indy,indx)=currs_plot(kk,3);
    end
    
    subplot(2,2,comp);set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    if (fl_plot_surf_imagesc == 1)
        surf(XX,YY,curr_on_grid);colormap('jet');colorbar; shading interp
        hold on;
        %plot3(currs_plot(:,1),currs_plot(:,2),abs(currs_plot(:,3)),'.','MarkerSize',15);
        view(0,90); axis tight; 
        if (slct_plane == 'xy');xlabel('x'); ylabel('y');
        elseif (slct_plane == 'xz');xlabel('x'); ylabel('z');
        elseif (slct_plane == 'yz');xlabel('y'); ylabel('z');
        end
        set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    elseif (fl_plot_surf_imagesc == 2)
        
        if (fl_ind_or_phys_len==1);
            x_var=XX(1,:);y_var=YY(:,1);
        else
            x_var=[1:1:length(XX(1,:))]; y_var=[1:1:length(YY(:,1))];
        end
        
        imagesc(x_var,y_var,curr_on_grid); colorbar; 
        if (slct_plane == 'xy');xlabel('x'); ylabel('y');
        elseif (slct_plane == 'xz');xlabel('x'); ylabel('z');
        elseif (slct_plane == 'yz');xlabel('y'); ylabel('z');
        end
        set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    end
    
    if (comp==1);title(str_Jx);
    elseif(comp==2);title(str_Jy);
    elseif(comp==3);title(str_Jz);
    end
    
end

