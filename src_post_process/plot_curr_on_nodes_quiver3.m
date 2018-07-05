function plot_curr_on_nodes_quiver3(nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned)

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

tmp_zeros_x=zeros(size(nodes_w_currs_x_aligned,1),1);
tmp_zeros_y=zeros(size(nodes_w_currs_y_aligned,1),1);
tmp_zeros_z=zeros(size(nodes_w_currs_z_aligned,1),1);

figure
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
q=quiver3([nodes_w_currs_x_aligned(:,1); nodes_w_currs_y_aligned(:,1); nodes_w_currs_z_aligned(:,1);],...
    [nodes_w_currs_x_aligned(:,2);nodes_w_currs_y_aligned(:,2);nodes_w_currs_z_aligned(:,2);],...
    [nodes_w_currs_x_aligned(:,3);nodes_w_currs_y_aligned(:,3);nodes_w_currs_z_aligned(:,3);],...
    [real(nodes_w_currs_x_aligned(:,4)); tmp_zeros_y;tmp_zeros_z;],...
    [tmp_zeros_x; real(nodes_w_currs_y_aligned(:,4));tmp_zeros_z;],...
    [tmp_zeros_x; tmp_zeros_y; real(nodes_w_currs_z_aligned(:,4));]);
xlabel('x'); ylabel('y');zlabel('z');
axis tight
set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');

q.Color = 'red';
q.MaxHeadSize = 5.0;
q.AutoScale = 'on';
q.AutoScaleFactor=15.0;
