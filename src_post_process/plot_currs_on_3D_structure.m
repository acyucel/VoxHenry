function plot_currs_on_3D_structure(x,Ae_only_leaving,r,Mc,dx)

interp_on=0; % use interpolated currents by mapping from centers to vertices
log_on = 0 ;% computed dB after normalizing with maximum current
dB_down = 50;
cut_on = 0; % take a cut on the geometry
cut_axis=3; %1->x,2->y,3->z, cut parallel to axis x,y,or z
cut_pnt=1e-6; % pnt through which cut passes
cut_below=1; % visualize below the cut [1] (<cut_pnt) or above the cut [0] (>cut_pnt)
fl_no_edge_line=1; % visualize without voxel edge lines
fl_title_on = 0; % add title to the figure

% 1) Obtain leaving currents on the nodes
num_curr=size(Ae_only_leaving,2);
x_node=Ae_only_leaving*x(1:num_curr);

%2) Retrieve voxel centers on which current exists
[L,M,N] = size(Mc); % domain size
tola=1e-12;
nonair_voxel_cens=zeros(3*num_curr/5,3);
dum=1;
for mz = 1:N;
    for my = 1:M;
        for mx = 1:L;
            if(abs(Mc(mx,my,mz)) > tola)
                nonair_voxel_cens(dum,1:3)=squeeze(r(mx,my,mz,1:3));
                nonair_voxel_cens(num_curr/5+dum,1:3)=squeeze(r(mx,my,mz,1:3));
                nonair_voxel_cens(2*num_curr/5+dum,1:3)=squeeze(r(mx,my,mz,1:3));
                dum=dum+1;
            end
        end;
    end;
end

%3) Find the currents on nodes - we are just interested in leaving currents

Ae_finding_nodes=Ae_only_leaving;
Ae_finding_nodes(:,num_curr/5*3+1:num_curr)=0;
[rows,cols,~] = find(Ae_finding_nodes);

%4) currents with voxels using exiting currents from nodes
voxels_w_currs=zeros(num_curr/5,4);
voxel_id=zeros(length(rows),1);
for kk=1:length(rows)
    %vox_cen_tmp=nonair_voxel_cens(cols(kk),1:3);
    if (cols(kk) > 0 && cols(kk) <= num_curr/5) % Jx current
        voxel_id(kk)=cols(kk);
    elseif (cols(kk) > num_curr/5 && cols(kk) <= 2*num_curr/5) % Jy current
        voxel_id(kk)=cols(kk)-(num_curr/5);
    elseif (cols(kk) > 2*num_curr/5 && cols(kk) <= 3*num_curr/5) % Jz current
        voxel_id(kk)=cols(kk)-(2*num_curr/5);
    end
    %voxels_w_currs(voxel_id,1:3)=nonair_voxel_cens(cols(kk),1:3);
    %voxels_w_currs(voxel_id,4)=voxels_w_currs(voxel_id,4)+x_node(rows(kk))^2;
end

voxels_w_currs(voxel_id(1:length(rows)),1:3)=nonair_voxel_cens(cols(1:length(rows)),1:3);

for kk=1:length(rows)
    voxels_w_currs(voxel_id(kk),4)=voxels_w_currs(voxel_id(kk),4)+x_node(rows(kk))^2;
end
%voxels_w_currs(voxel_id(1:length(rows)),4)=voxels_w_currs(voxel_id(1:length(rows)),4)+x_node(rows((1:length(rows)))).^2;

% since we visualize magnitude of current:
voxels_w_currs(:,4)=abs(sqrt(voxels_w_currs(:,4)));

% creating faces of voxels and assigning value to them

faces_w_currs_x=zeros(4,6*num_curr/5); % each voxel faces - vertices x coordinate
faces_w_currs_y=zeros(4,6*num_curr/5); % each voxel faces - vertices y coordinate
faces_w_currs_z=zeros(4,6*num_curr/5); % each voxel faces - vertices z coordinate
faces_w_currs_color=zeros(1,6*num_curr/5); % each voxel faces' color

half_dx=dx/2;
tola=0; %1e-12;
num_vox_w_currs=size(voxels_w_currs,1);
for kk=1:num_vox_w_currs
    color_tmp=voxels_w_currs(kk,4);
    vox_cen_tmp=voxels_w_currs(kk,1:3);
    % left face
    faces_w_currs_x(1:4,kk)=vox_cen_tmp(1)-half_dx+tola;
    faces_w_currs_y(1:2,kk)=vox_cen_tmp(2)-half_dx;
    faces_w_currs_y(3:4,kk)=vox_cen_tmp(2)+half_dx;
    faces_w_currs_z([2 3],kk)=vox_cen_tmp(3)-half_dx;
    faces_w_currs_z([1 4],kk)=vox_cen_tmp(3)+half_dx;
    
    % right face
    ind_dum=num_vox_w_currs+kk;
    faces_w_currs_x(1:4,ind_dum)=vox_cen_tmp(1)+half_dx-tola;
    faces_w_currs_y(1:2,ind_dum)=vox_cen_tmp(2)-half_dx;
    faces_w_currs_y(3:4,ind_dum)=vox_cen_tmp(2)+half_dx;
    faces_w_currs_z([2 3],ind_dum)=vox_cen_tmp(3)-half_dx;
    faces_w_currs_z([1 4],ind_dum)=vox_cen_tmp(3)+half_dx;
    
    % front face
    ind_dum=2*num_vox_w_currs+kk;
    faces_w_currs_x(1:2,ind_dum)=vox_cen_tmp(1)-half_dx;
    faces_w_currs_x(3:4,ind_dum)=vox_cen_tmp(1)+half_dx;
    
    faces_w_currs_y(1:4,ind_dum)=vox_cen_tmp(2)-half_dx+tola;
    
    faces_w_currs_z([2 3],ind_dum)=vox_cen_tmp(3)-half_dx;
    faces_w_currs_z([1 4],ind_dum)=vox_cen_tmp(3)+half_dx;
    
    % back face
    ind_dum=3*num_vox_w_currs+kk;
    faces_w_currs_x(1:2,ind_dum)=vox_cen_tmp(1)-half_dx;
    faces_w_currs_x(3:4,ind_dum)=vox_cen_tmp(1)+half_dx;
    
    faces_w_currs_y(1:4,ind_dum)=vox_cen_tmp(2)+half_dx-tola;
    
    faces_w_currs_z([2 3],ind_dum)=vox_cen_tmp(3)-half_dx;
    faces_w_currs_z([1 4],ind_dum)=vox_cen_tmp(3)+half_dx;
    
    % bottom face
    ind_dum=4*num_vox_w_currs+kk;
    faces_w_currs_x(1:2,ind_dum)=vox_cen_tmp(1)-half_dx;
    faces_w_currs_x(3:4,ind_dum)=vox_cen_tmp(1)+half_dx;
    
    faces_w_currs_y([2 3],ind_dum)=vox_cen_tmp(2)-half_dx;
    faces_w_currs_y([1 4],ind_dum)=vox_cen_tmp(2)+half_dx;
    
    faces_w_currs_z(1:4,ind_dum)=vox_cen_tmp(3)-half_dx;
    
    % top face
    ind_dum=5*num_vox_w_currs+kk;
    
    faces_w_currs_x(1:2,ind_dum)=vox_cen_tmp(1)-half_dx;
    faces_w_currs_x(3:4,ind_dum)=vox_cen_tmp(1)+half_dx;
    
    faces_w_currs_y([2 3],ind_dum)=vox_cen_tmp(2)-half_dx;
    faces_w_currs_y([1 4],ind_dum)=vox_cen_tmp(2)+half_dx;
    
    faces_w_currs_z(1:4,ind_dum)=vox_cen_tmp(3)+half_dx;
    
    % assign colors
    faces_w_currs_color(1,[0 num_vox_w_currs 2*num_vox_w_currs ...
        3*num_vox_w_currs 4*num_vox_w_currs 5*num_vox_w_currs]+kk)=color_tmp;
    
end

if (interp_on == 1) % we compute the currents on vertices by averaging
    % interpolated visualization
    
    max_ind_x=max(max(faces_w_currs_x(1:4,:)/dx));
    min_ind_x=min(min(faces_w_currs_x(1:4,:)/dx));
    
    max_ind_y=max(max(faces_w_currs_y(1:4,:)/dx));
    min_ind_y=min(min(faces_w_currs_y(1:4,:)/dx));
    
    max_ind_z=max(max(faces_w_currs_z(1:4,:)/dx));
    min_ind_z=min(min(faces_w_currs_z(1:4,:)/dx));
    
    loc_org=[min_ind_x min_ind_y min_ind_z]*dx-[dx dx dx];
    
    tmp_tensor=zeros(round((max_ind_x-min_ind_x))+1,round((max_ind_y-min_ind_y))+1,round((max_ind_z-min_ind_z))+1,2);
    for kk=1:size(faces_w_currs_x,2)
        for ll=1:4
            tmp_coor=[faces_w_currs_x(ll,kk) faces_w_currs_y(ll,kk) faces_w_currs_z(ll,kk)];
            ind_xyz=round((tmp_coor-loc_org)/dx);
            tmp_tensor(ind_xyz(1),ind_xyz(2),ind_xyz(3),1)=tmp_tensor(ind_xyz(1),ind_xyz(2),ind_xyz(3),1)+...
                faces_w_currs_color(kk);
            tmp_tensor(ind_xyz(1),ind_xyz(2),ind_xyz(3),2)=tmp_tensor(ind_xyz(1),ind_xyz(2),ind_xyz(3),2)+1;
        end
    end
    
    for kk=1:size(tmp_tensor,1)
        for ll=1:size(tmp_tensor,2)
            for mm=1:size(tmp_tensor,3)
                if (tmp_tensor(kk,ll,mm,1) ~= 0)
                    tmp_tensor(kk,ll,mm,1)=tmp_tensor(kk,ll,mm,1)/tmp_tensor(kk,ll,mm,2);
                end
            end
        end
    end
    
    % mapping back to vector from tensor
    
    vertices_w_currs_color=zeros(4,size(faces_w_currs_x,2));
    
    for kk=1:size(faces_w_currs_x,2)
        for ll=1:4
            tmp_coor=[faces_w_currs_x(ll,kk) faces_w_currs_y(ll,kk) faces_w_currs_z(ll,kk)];
            ind_xyz=round((tmp_coor-loc_org)/dx);
            vertices_w_currs_color(ll,kk)=tmp_tensor(ind_xyz(1),ind_xyz(2),ind_xyz(3),1);
        end
    end
    clear tmp_tensor
end


if (log_on == 1) % for normalization, finding the maximum current
    if (interp_on == 1)
        maxx_curr=max(max(vertices_w_currs_color));
    else
        maxx_curr=max(max(faces_w_currs_color));
    end
    
end

if (cut_on == 1)
    if(cut_axis == 1) % cut parallel to x
        faces_w_currs_tmp=faces_w_currs_x;
    elseif (cut_axis == 2) % cut parallel to y
        faces_w_currs_tmp=faces_w_currs_y;
    elseif (cut_axis == 3) % cut parallel to z
        faces_w_currs_tmp=faces_w_currs_z;
    end
    if (cut_below==1) % visualize below the cut
        counter=0;
        for kk=1:size(faces_w_currs_tmp,2)
            if (faces_w_currs_tmp(1,kk) <= cut_pnt && faces_w_currs_tmp(2,kk) <= cut_pnt && ...
                    faces_w_currs_tmp(3,kk) <= cut_pnt && faces_w_currs_tmp(4,kk) <= cut_pnt)
                counter=counter+1;
            end
        end
        
        inds_cut=zeros(counter,1);
        
        dum=1;
        for kk=1:size(faces_w_currs_tmp,2)
            if (faces_w_currs_tmp(1,kk) <= cut_pnt && faces_w_currs_tmp(2,kk) <= cut_pnt && ...
                    faces_w_currs_tmp(3,kk) <= cut_pnt && faces_w_currs_tmp(4,kk) <= cut_pnt)
                inds_cut(dum) = kk;
                dum=dum+1;
            end
        end
    else % visualize above the cut
        counter=0;
        for kk=1:size(faces_w_currs_tmp,2)
            if (faces_w_currs_tmp(1,kk) >= cut_pnt && faces_w_currs_tmp(2,kk) >= cut_pnt && ...
                    faces_w_currs_tmp(3,kk) >= cut_pnt && faces_w_currs_tmp(4,kk) >= cut_pnt)
                counter=counter+1;
            end
        end
        
        inds_cut=zeros(counter,1);
        
        dum=1;
        for kk=1:size(faces_w_currs_tmp,2)
            if (faces_w_currs_tmp(1,kk) >= cut_pnt && faces_w_currs_tmp(2,kk) >= cut_pnt && ...
                    faces_w_currs_tmp(3,kk) >= cut_pnt && faces_w_currs_tmp(4,kk) >= cut_pnt)
                inds_cut(dum) = kk;
                dum=dum+1;
            end
        end
    end
    
    figure
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    if(interp_on == 1)
        if (log_on ==1)
            tmp_vect=[];
            tmp_vect=max(20*log10(vertices_w_currs_color(:,inds_cut)/maxx_curr),-dB_down);
            ppp=patch(faces_w_currs_x(:,inds_cut),faces_w_currs_y(:,inds_cut),faces_w_currs_z(:,inds_cut),tmp_vect);
            colorbar; caxis([-dB_down 0]);
            if(fl_title_on==1); title('Normalized Total Current (dB)'); end
        else
            ppp=patch(faces_w_currs_x(:,inds_cut),faces_w_currs_y(:,inds_cut),faces_w_currs_z(:,inds_cut),vertices_w_currs_color(:,inds_cut));
            colorbar;
            if(fl_title_on==1); title('Total Current'); end;
        end
    else
        if (log_on ==1)
            tmp_vect=[];
            tmp_vect=max(20*log10(faces_w_currs_color(inds_cut)/maxx_curr),-dB_down);
            ppp=patch(faces_w_currs_x(:,inds_cut),faces_w_currs_y(:,inds_cut),faces_w_currs_z(:,inds_cut),tmp_vect);
            colorbar; caxis([-dB_down 0]);
            if(fl_title_on==1); title('Normalized Total Current (dB)'); end;
        else
            ppp=patch(faces_w_currs_x(:,inds_cut),faces_w_currs_y(:,inds_cut),faces_w_currs_z(:,inds_cut),faces_w_currs_color(inds_cut));
            colorbar
            if(fl_title_on==1); title('Total Current'); end;
        end
    end
    map='jet';
    colormap(map);
    axis tight; grid on; % this should be here, otherwise it doesn't add vertical patches
    axis equal; grid on;
    view(130,20)
    
    xlabel('x');ylabel('y');zlabel('z');
    set(gcf,'renderer','openGL');
    
else
    FigHandle = figure;
    %set(FigHandle, 'Position', [100, 100, 1280, 1024]);
    set(FigHandle, 'Position', [100, 100, 1280, 500]);
    %subplot(2,1,1)
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    if(interp_on == 1)
        if (log_on ==1)
            tmp_vect=[];
            tmp_vect=max(20*log10(vertices_w_currs_color/maxx_curr),-dB_down);
            ppp=patch(faces_w_currs_x,faces_w_currs_y,faces_w_currs_z,tmp_vect);
            colorbar; caxis([-dB_down 0]);
            if(fl_title_on==1); title('Normalized Total Current (dB)'); end;
        else
            ppp=patch(faces_w_currs_x,faces_w_currs_y,faces_w_currs_z,vertices_w_currs_color);
            colorbar
            if(fl_title_on==1); title('Total Current'); end
        end
    else
        if (log_on ==1)
            tmp_vect=[];
            tmp_vect=max(20*log10(faces_w_currs_color/maxx_curr),-dB_down);
            ppp=patch(faces_w_currs_x,faces_w_currs_y,faces_w_currs_z,tmp_vect);
            colorbar; caxis([-dB_down 0]);
            if(fl_title_on==1); title('Normalized Total Current (dB)'); end;
        else
            ppp=patch(faces_w_currs_x,faces_w_currs_y,faces_w_currs_z,faces_w_currs_color);
            colorbar
            if(fl_title_on==1); title('Total Current'); end;
        end
    end
    map='hot';
    %map='bone';
    hcb=colormap(hot);%caxis([0 20.0e-3]); 
    %set(hcb,'XTick',[0 5e-3,10e-3,18e-3])
    axis tight; grid on; % this should be here, otherwise it doesn't add vertical patches
    axis equal; grid on;
    %view(130,20)
    view(-45,20)
    xlabel('x');ylabel('y');zlabel('z');
    %set(gcf,'Renderer','opengl');
    set(gcf,'Renderer','zbuffer')
    set(gca,'FontSize',24); set(gca,'FontName','Times New Roman');
    %colorbar off
    colorbar
    
end

%set(ppp,'facealpha',1.0);
%set(ppp,'FaceColor','interp');
if (fl_no_edge_line == 1)
    set(ppp,'EdgeColor','none')
end