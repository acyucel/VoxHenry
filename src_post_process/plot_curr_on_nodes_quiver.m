function plot_curr_on_nodes_quiver(slct_plane,slct_cut,nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned)

fl_slct_part=1; %1-real,2-imag,3-abs

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

tola=1e-12;
if (slct_plane == 'xy') % we have x and y aligned nodes
    tmp_book_keep_x=zeros(size(nodes_w_currs_x_aligned,1),1);
    tmp_book_keep_y=zeros(size(nodes_w_currs_y_aligned,1),1);
    
    % find the x aligned nodes on specified cut
    for kk=1:size(nodes_w_currs_x_aligned,1)
        if (abs(nodes_w_currs_x_aligned(kk,3)-slct_cut)<tola)
            tmp_book_keep_x(kk)=1;
        end
    end
    
    % find the y aligned nodes on specified cut
    for kk=1:size(nodes_w_currs_y_aligned,1)
        if (abs(nodes_w_currs_y_aligned(kk,3)-slct_cut)<tola)
            tmp_book_keep_y(kk)=1;
        end
    end
    
    % put the nodes in seperate array for plotting
    currs_x_plot= zeros(sum(tmp_book_keep_x),4);
    dum=1;
    for kk=1:size(nodes_w_currs_x_aligned,1)
        if (tmp_book_keep_x(kk) == 1)
            currs_x_plot(dum,1:2)=nodes_w_currs_x_aligned(kk,1:2);
            currs_x_plot(dum,3)=nodes_w_currs_x_aligned(kk,4);
            dum=dum+1;
        end
    end
    
    currs_y_plot= zeros(sum(tmp_book_keep_y),4);
    dum=1;
    for kk=1:size(nodes_w_currs_y_aligned,1)
        if (tmp_book_keep_y(kk) == 1)
            currs_y_plot(dum,1:2)=nodes_w_currs_y_aligned(kk,1:2);
            currs_y_plot(dum,4)=nodes_w_currs_y_aligned(kk,4);
            dum=dum+1;
        end
    end
    
    figure;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    q = quiver([currs_x_plot(:,1);currs_y_plot(:,1)],...
        [currs_x_plot(:,2);currs_y_plot(:,2)],...
        [currs_x_plot(:,3);currs_y_plot(:,3)],...
        [currs_x_plot(:,4);currs_y_plot(:,4)]);
    xlabel('x'); ylabel('y');title([str_Jx,' & ',str_Jy]);
    axis tight
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
end

tola=1e-12;
if (slct_plane == 'yz') % we have y and z aligned nodes
    tmp_book_keep_y=zeros(size(nodes_w_currs_y_aligned,1),1);
    tmp_book_keep_z=zeros(size(nodes_w_currs_z_aligned,1),1);
    
    % find the z aligned nodes on specified cut
    for kk=1:size(nodes_w_currs_z_aligned,1)
        if (abs(nodes_w_currs_z_aligned(kk,1)-slct_cut)<tola)
            tmp_book_keep_z(kk)=1;
        end
    end
    
    % find the y aligned nodes on specified cut
    for kk=1:size(nodes_w_currs_y_aligned,1)
        if (abs(nodes_w_currs_y_aligned(kk,1)-slct_cut)<tola)
            tmp_book_keep_y(kk)=1;
        end
    end
    
    % put the nodes in seperate array for plotting
    currs_z_plot= zeros(sum(tmp_book_keep_z),4);
    dum=1;
    for kk=1:size(nodes_w_currs_z_aligned,1)
        if (tmp_book_keep_z(kk) == 1)
            % first entry y, second entry z
            currs_z_plot(dum,1:2)=nodes_w_currs_z_aligned(kk,2:3); 
            currs_z_plot(dum,4)=nodes_w_currs_z_aligned(kk,4);
            dum=dum+1;
        end
    end
    
    currs_y_plot= zeros(sum(tmp_book_keep_y),4);
    dum=1;
    for kk=1:size(nodes_w_currs_y_aligned,1)
        if (tmp_book_keep_y(kk) == 1)
            currs_y_plot(dum,1:2)=nodes_w_currs_y_aligned(kk,2:3);
            currs_y_plot(dum,3)=nodes_w_currs_y_aligned(kk,4);
            dum=dum+1;
        end
    end
    
    figure;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    q = quiver([currs_y_plot(:,1);currs_z_plot(:,1)],...
        [currs_y_plot(:,2);currs_z_plot(:,2)],...
        [currs_y_plot(:,3);currs_z_plot(:,3)],...
        [currs_y_plot(:,4);currs_z_plot(:,4)]);
    xlabel('y'); ylabel('z'); title([str_Jy,' & ',str_Jz]);
    axis tight
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    
end

tola=1e-12;
if (slct_plane == 'xz') % we have y and z aligned nodes
    tmp_book_keep_x=zeros(size(nodes_w_currs_x_aligned,1),1);
    tmp_book_keep_z=zeros(size(nodes_w_currs_z_aligned,1),1);
    
    % find the z aligned nodes on specified cut
    for kk=1:size(nodes_w_currs_z_aligned,1)
        if (abs(nodes_w_currs_z_aligned(kk,2)-slct_cut)<tola)
            tmp_book_keep_z(kk)=1;
        end
    end
    
    % find the x aligned nodes on specified cut
    for kk=1:size(nodes_w_currs_x_aligned,1)
        if (abs(nodes_w_currs_x_aligned(kk,2)-slct_cut)<tola)
            tmp_book_keep_x(kk)=1;
        end
    end
    
    % put the nodes in seperate array for plotting
    currs_z_plot= zeros(sum(tmp_book_keep_z),4);
    dum=1;
    for kk=1:size(nodes_w_currs_z_aligned,1)
        if (tmp_book_keep_z(kk) == 1)
            % first entry x, second entry z
            currs_z_plot(dum,1:2)=nodes_w_currs_z_aligned(kk,[1 3]); 
            currs_z_plot(dum,4)=nodes_w_currs_z_aligned(kk,4);
            dum=dum+1;
        end
    end
    
    currs_x_plot= zeros(sum(tmp_book_keep_x),4);
    dum=1;
    for kk=1:size(nodes_w_currs_x_aligned,1)
        if (tmp_book_keep_x(kk) == 1)
            currs_x_plot(dum,1:2)=nodes_w_currs_x_aligned(kk,[1 3]);
            currs_x_plot(dum,3)=nodes_w_currs_x_aligned(kk,4);
            dum=dum+1;
        end
    end
    
    figure;
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    q = quiver([currs_x_plot(:,1);currs_z_plot(:,1)],...
        [currs_x_plot(:,2);currs_z_plot(:,2)],...
        [currs_x_plot(:,3);currs_z_plot(:,3)],...
        [currs_x_plot(:,4);currs_z_plot(:,4)]);
    xlabel('x'); ylabel('z'); title([str_Jx,' & ',str_Jz]);
    axis tight
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    
end

q.Color = 'red';
q.MaxHeadSize = 5.0;
q.AutoScale = 'on';
q.AutoScaleFactor=5.0;
