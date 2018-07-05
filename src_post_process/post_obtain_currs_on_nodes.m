function [nodes_w_currs_x_aligned,nodes_w_currs_y_aligned,nodes_w_currs_z_aligned]=post_obtain_currs_on_nodes(x,Ae_only_leaving,Ae_only_entering_bndry,r,Mc,dx)

fl_cpu_profile=0;
num_curr=size(Ae_only_leaving,2);
% 1) Obtain currents on the nodes
if (fl_cpu_profile==1); tic; end;
x_inner_node=Ae_only_leaving*x(1:num_curr);
x_bndry_node=Ae_only_entering_bndry*x(1:num_curr);
%x_node=x_inner_node+x_bndry_node;
x_node=-x_inner_node+x_bndry_node;
if (fl_cpu_profile==1); disp(['Time for obtaining currents on nodes ::: ',num2str(toc)]); end;

%2) Retrieve voxel centers on which current exists
if (fl_cpu_profile==1); tic; end;
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
if (fl_cpu_profile==1); disp(['Time for retrieving voxel centers ::: ',num2str(toc)]); end;

%3) Find the currents on nodes

if (fl_cpu_profile==1); tic; end;
Ae_finding_nodes=Ae_only_leaving;
Ae_finding_nodes(:,num_curr/5*3+1:num_curr)=0;
Ae_finding_nodes2=Ae_only_entering_bndry;
Ae_finding_nodes2(:,num_curr/5*3+1:num_curr)=0;

[rows,cols,vals] = find(Ae_finding_nodes);
[rows2,cols2,vals2] = find(Ae_finding_nodes2);
 
if (fl_cpu_profile==1); disp(['Time for getting modified Aes::: ',num2str(toc)]); end;

if (fl_cpu_profile==1); tic; end;
nodes_w_currs=zeros(length(rows)+length(rows2),5);

% for exiting currents
for kk=1:length(rows)
     vox_cen_tmp=nonair_voxel_cens(cols(kk),1:3);
     if (cols(kk) > 0 && cols(kk) <= num_curr/5) % x aligned panel
         %nodes_w_currs(kk,1:3)=[vox_cen_tmp(1)-dx/2 vox_cen_tmp(2) vox_cen_tmp(3)];
         vox_cen_tmp(1)=vox_cen_tmp(1)-dx/2;
     elseif (cols(kk) > num_curr/5 && cols(kk) <= 2*num_curr/5) % y aligned panel
         %nodes_w_currs(kk,1:3)=[vox_cen_tmp(1) vox_cen_tmp(2)-dx/2 vox_cen_tmp(3)];
         vox_cen_tmp(2)=vox_cen_tmp(2)-dx/2;
     elseif (cols(kk) > 2*num_curr/5 && cols(kk) <= 3*num_curr/5) % z aligned panel
         vox_cen_tmp(3)=vox_cen_tmp(3)-dx/2;
         %nodes_w_currs(kk,1:3)=[vox_cen_tmp(1) vox_cen_tmp(2) vox_cen_tmp(3)-dx/2];
     end
     nodes_w_currs(kk,1:3) = vox_cen_tmp;
     %nodes_w_currs(kk,4) = x_node(rows(kk));
     %nodes_w_currs(kk,5) = rows(kk); % original node id
end
% Attention: the following is moved here for accelerating do-loop
nodes_w_currs(1:length(rows),4) = x_node(rows(1:length(rows)));
nodes_w_currs(1:length(rows),5) = rows(1:length(rows)); % original node id

if (fl_cpu_profile==1); disp(['Time for nodes w/ exiting currs ::: ',num2str(toc)]); end;

% for entering currents to boundaries
if (fl_cpu_profile==1); tic; end;
tmp_vect=zeros(length(rows2),3);
for kk=1:length(rows2)
    vox_cen_tmp=nonair_voxel_cens(cols2(kk),1:3);
    if (cols2(kk) > 0 && cols2(kk) <= num_curr/5) % x aligned panel
        %nodes_w_currs(length(rows)+kk,1:3)=[vox_cen_tmp(1)+dx/2 vox_cen_tmp(2) vox_cen_tmp(3)];
        vox_cen_tmp(1)=vox_cen_tmp(1)+dx/2;
    elseif (cols2(kk) > num_curr/5 && cols2(kk) <= 2*num_curr/5) % y aligned panel
        %nodes_w_currs(length(rows)+kk,1:3)=[vox_cen_tmp(1) vox_cen_tmp(2)+dx/2 vox_cen_tmp(3)];
        vox_cen_tmp(2)=vox_cen_tmp(2)+dx/2;
    elseif (cols2(kk) > 2*num_curr/5 && cols2(kk) <= 3*num_curr/5) % z aligned panel
        %nodes_w_currs(length(rows)+kk,1:3)=[vox_cen_tmp(1) vox_cen_tmp(2) vox_cen_tmp(3)+dx/2];
        vox_cen_tmp(3)=vox_cen_tmp(3)+dx/2;
    end
    tmp_vect(kk,1:3)=vox_cen_tmp;
    %nodes_w_currs(length(rows)+kk,1:3)=vox_cen_tmp;
    %nodes_w_currs(length(rows)+kk,4)=x_node(rows2(kk));
    %nodes_w_currs(length(rows)+kk,5)=rows2(kk);
end
% Attention: the following is moved here for accelerating do-loop
%ind_dum=[length(rows)+1:1:length(rows)+length(rows2)];
ind_dum=length(rows)+1:1:length(rows)+length(rows2);
nodes_w_currs(ind_dum,1:3)=tmp_vect;
nodes_w_currs(ind_dum,4)=x_node(rows2(1:length(rows2)));
nodes_w_currs(ind_dum,5)=rows2(1:length(rows2));
clear ind_dum tmp_vect

if (fl_cpu_profile==1); disp(['Time for nodes w/ entering currs ::: ',num2str(toc)]); end;

% sorting nodes
if (fl_cpu_profile==1); tic; end;
nodes_w_currs_sorted=zeros(size(nodes_w_currs,1),4);

% for kk=1:size(nodes_w_currs,1)
%     nodes_w_currs_sorted(round(nodes_w_currs(kk,5)),1:4)=nodes_w_currs(kk,1:4);
% end
% Attention: the following is moved here for accelerating do-loop
nodes_w_currs_sorted(round(nodes_w_currs(1:size(nodes_w_currs,1),5)),1:4)=nodes_w_currs(1:size(nodes_w_currs,1),1:4);


if (fl_cpu_profile==1); disp(['Time for sorting currs and nodes ::: ',num2str(toc)]); end;

clear nodes_w_currs

% 4) Classify the nodes & currents for defining currents w/ components
% x-aligned nodes has Jx, y-aligned nodes Jy, and z-aligned nodes Jz
% x-aligned nodes
if (fl_cpu_profile==1); tic; end;
tmp_vect=zeros(num_curr,1);tmp_vect(1:num_curr/5)=1;
x_aligned_nodes_exit=find(Ae_finding_nodes*tmp_vect);
x_aligned_nodes_enter=find(Ae_finding_nodes2*tmp_vect);
if (fl_cpu_profile==1); disp(['Time for getting x-aligned nodes w/currs::: ',num2str(toc)]); end;
% y-aligned nodes
if (fl_cpu_profile==1); tic; end;
tmp_vect=zeros(num_curr,1);tmp_vect(num_curr/5+1:2*num_curr/5)=1;
y_aligned_nodes_exit=find(Ae_finding_nodes*tmp_vect);
y_aligned_nodes_enter=find(Ae_finding_nodes2*tmp_vect);
if (fl_cpu_profile==1); disp(['Time for getting y-aligned nodes w/currs::: ',num2str(toc)]); end;
% z-aligned nodes
if (fl_cpu_profile==1); tic; end;
tmp_vect=zeros(num_curr,1);tmp_vect(2*num_curr/5+1:3*num_curr/5)=1;
z_aligned_nodes_exit=find(Ae_finding_nodes*tmp_vect);
z_aligned_nodes_enter=find(Ae_finding_nodes2*tmp_vect);
if (fl_cpu_profile==1); disp(['Time for getting z-aligned nodes w/currs::: ',num2str(toc)]); end;
%clear Ae_finding_nodes Ae_finding_nodes2

tot_num_x_nodes=length(x_aligned_nodes_exit)+length(x_aligned_nodes_enter);
tot_num_y_nodes=length(y_aligned_nodes_exit)+length(y_aligned_nodes_enter);
tot_num_z_nodes=length(z_aligned_nodes_exit)+length(z_aligned_nodes_enter);

nodes_w_currs_x_aligned=zeros(tot_num_x_nodes,4);
nodes_w_currs_y_aligned=zeros(tot_num_y_nodes,4);
nodes_w_currs_z_aligned=zeros(tot_num_z_nodes,4);

% x-aligned nodes with currents
if (fl_cpu_profile==1); tic; end;
%for kk=1:length(x_aligned_nodes_exit)
%    nodes_w_currs_x_aligned(kk,:) = nodes_w_currs_sorted(x_aligned_nodes_exit(kk),:);
%end
% Attention: the following is moved here for accelerating do-loop
nodes_w_currs_x_aligned(1:length(x_aligned_nodes_exit),:) = nodes_w_currs_sorted(x_aligned_nodes_exit(1:length(x_aligned_nodes_exit)),:);

%for kk=length(x_aligned_nodes_exit)+1:tot_num_x_nodes
%    nodes_w_currs_x_aligned(kk,:) = nodes_w_currs_sorted(x_aligned_nodes_enter(kk-length(x_aligned_nodes_exit)),:);
%end
% Attention: the following is moved here for accelerating do-loop
nodes_w_currs_x_aligned([length(x_aligned_nodes_exit)+1:tot_num_x_nodes],:) = nodes_w_currs_sorted(x_aligned_nodes_enter([length(x_aligned_nodes_exit)+1:tot_num_x_nodes]-length(x_aligned_nodes_exit)),:);


if (fl_cpu_profile==1); disp(['Time for obtaining nodes_w_currs_x_aligned ::: ',num2str(toc)]); end;

% y-aligned nodes with currents
if (fl_cpu_profile==1); tic; end;
% for kk=1:length(y_aligned_nodes_exit)
%     nodes_w_currs_y_aligned(kk,:) = nodes_w_currs_sorted(y_aligned_nodes_exit(kk),:);
% end
% Attention: the following is moved here for accelerating do-loop
nodes_w_currs_y_aligned(1:length(y_aligned_nodes_exit),:) = nodes_w_currs_sorted(y_aligned_nodes_exit(1:length(y_aligned_nodes_exit)),:);

%for kk=length(y_aligned_nodes_exit)+1:tot_num_y_nodes
%    nodes_w_currs_y_aligned(kk,:) = nodes_w_currs_sorted(y_aligned_nodes_enter(kk-length(y_aligned_nodes_exit)),:);
%end
% Attention: the following is moved here for accelerating do-loop
nodes_w_currs_y_aligned([length(y_aligned_nodes_exit)+1:tot_num_y_nodes],:) = nodes_w_currs_sorted(y_aligned_nodes_enter([length(y_aligned_nodes_exit)+1:tot_num_y_nodes]-length(y_aligned_nodes_exit)),:);

if (fl_cpu_profile==1); disp(['Time for obtaining nodes_w_currs_y_aligned ::: ',num2str(toc)]); end;

% z-aligned nodes with currents
if (fl_cpu_profile==1); tic; end;
%for kk=1:length(z_aligned_nodes_exit)
%    nodes_w_currs_z_aligned(kk,:) = nodes_w_currs_sorted(z_aligned_nodes_exit(kk),:);
%end
% Attention: the following is moved here for accelerating do-loop
nodes_w_currs_z_aligned(1:length(z_aligned_nodes_exit),:) = nodes_w_currs_sorted(z_aligned_nodes_exit(1:length(z_aligned_nodes_exit)),:);

%for kk=length(z_aligned_nodes_exit)+1:tot_num_z_nodes
%    nodes_w_currs_z_aligned(kk,:) = nodes_w_currs_sorted(z_aligned_nodes_enter(kk-length(z_aligned_nodes_exit)),:);
%end
% Attention: the following is moved here for accelerating do-loop
nodes_w_currs_z_aligned([length(z_aligned_nodes_exit)+1:tot_num_z_nodes],:) = ...
    nodes_w_currs_sorted(z_aligned_nodes_enter([length(z_aligned_nodes_exit)+1:tot_num_z_nodes]-length(z_aligned_nodes_exit)),:);

if (fl_cpu_profile==1); disp(['Time for obtaining nodes_w_currs_z_aligned ::: ',num2str(toc)]); end;

