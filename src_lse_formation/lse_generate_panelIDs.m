function [all_panels_ids, num_nodes] = lse_generate_panelIDs(idxS, L, M, N)

fl_profile = 0;

% constants
num_nonair_cube=length(idxS);


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
voxel2panel_x=zeros(num_nonair_cube,2); % 1 left, 2 right panel (i.e -x, +x)
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
voxel2panel_y=zeros(num_nonair_cube,2); % 1 front panel, 2 back panel (i.e. -y, +y)
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
voxel2panel_z=zeros(num_nonair_cube,2); % 1 bottom panel, 2 top panel (i.e. -z / +z)
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
% stitch the information together
%    

voxel2panel_y = num_panel_x + voxel2panel_y;

voxel2panel_z = (num_panel_x+num_panel_y) + voxel2panel_z;

% in a nutshell, 'all_panels_ids' is a matrix that has rows corresponding to the non-empty voxels,
% and columns giving the node IDs of the six nodes belonging to the voxel specified by that row,
% where nodes are in column order -x/+x/-y/+y/-z/+z 
all_panels_ids=[voxel2panel_x voxel2panel_y voxel2panel_z];
num_nodes=num_panel_x+num_panel_y+num_panel_z;

if(fl_profile == 1); disp(['Time for indexing panels ::: ', num2str(toc)]); end;


