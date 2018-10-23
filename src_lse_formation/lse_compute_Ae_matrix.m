function [Ae, Ae_only_leaving, Ae_only_entering_bndry] = lse_compute_Ae_matrix(idxS, all_panels_ids, num_nodes)

tstart = tic;
fl_profile = 0;

% constants
num_nonair_cube=length(idxS);

%
% Generating Ae matrix
% 

disp('  Generating Ae matrix by inspection...')


% forming Ae matrix

tic

% enter +1, leaves -1; first entry for leaving, second entry for entering
const_lin=1/2;

sp_mat_inds=zeros(16*num_nonair_cube,3);
sp_mat_inds_only_leaving_currs=zeros(8*num_nonair_cube,3);

% 'sp_mat_inds' contains elements in strict column order (even if the rows are not ordered,
% but they are sparse)

for kk=1:num_nonair_cube % pertinent to Jx, Jy, and Jz
    sp_mat_inds(kk,                  1:3) = [all_panels_ids(kk,1) kk -1];
    sp_mat_inds(kk+num_nonair_cube,  1:3) = [all_panels_ids(kk,2) kk 1];
    sp_mat_inds(kk+2*num_nonair_cube,1:3) = [all_panels_ids(kk,3) kk+num_nonair_cube -1];
    sp_mat_inds(kk+3*num_nonair_cube,1:3) = [all_panels_ids(kk,4) kk+num_nonair_cube 1];
    sp_mat_inds(kk+4*num_nonair_cube,1:3) = [all_panels_ids(kk,5) kk+2*num_nonair_cube -1];
    sp_mat_inds(kk+5*num_nonair_cube,1:3) = [all_panels_ids(kk,6) kk+2*num_nonair_cube 1];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(kk,1:3)= [all_panels_ids(kk,1) kk -1];
    sp_mat_inds_only_leaving_currs(kk+num_nonair_cube,1:3)= [all_panels_ids(kk,3) kk+num_nonair_cube -1];
    sp_mat_inds_only_leaving_currs(kk+2*num_nonair_cube,1:3)= [all_panels_ids(kk,5) kk+2*num_nonair_cube -1];
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
    sp_mat_inds(8*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,3) kk -const_lin];
    sp_mat_inds(9*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,4) kk -const_lin];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(4*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,3) kk -const_lin];
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
    sp_mat_inds(12*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,3) kk const_lin];
    sp_mat_inds(13*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,4) kk const_lin];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(6*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,3) kk const_lin];
    dum=dum+1;
end

dum=1;
for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsz
    sp_mat_inds(14*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,5) kk -2*const_lin];
    sp_mat_inds(15*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,6) kk -2*const_lin];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(7*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,5) kk -2*const_lin];
    dum=dum+1;
end

Ae=sparse(sp_mat_inds(:,1),sp_mat_inds(:,2),sp_mat_inds(:,3),num_nodes,5*num_nonair_cube);

infomem1 = whos('Ae');
memestimated = (infomem1.bytes)/(1024*1024);
disp(['  Memory for Ae matrix (MB)::' , num2str(memestimated)]);

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

disp('  Done... Generating Ae matrix by inspection...')

if(fl_profile == 1); disp(['Total time for generating Ae matrix ::: ', num2str(toc(tstart))]); end;


num_node = size(Ae,1);
num_curr = size(Ae,2);
disp('-----------------------------------------------------')
disp(['Number of current unknowns ::: ',num2str(num_curr)])
disp(['Number of potential unknowns ::: ',num2str(num_node)])
disp(['Number of total unknowns ::: ',num2str(num_curr+num_node)])
disp('-----------------------------------------------------')

