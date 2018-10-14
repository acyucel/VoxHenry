% -------------------------------------------------------------------------
%                 Define EM Vars/Constants and Domain Parameters
% -------------------------------------------------------------------------

pre_define_structure_params
dx = Res;

% ------------------------------------------------------------------------
%                  Obtain port Nodes
% -------------------------------------------------------------------------


% constants
num_nonair_cube=length(idxS);


% Get the voxel2nodes matrix (each row is a non-empty panel, each of the six columns
% is a node ID in the order -x/+x/-y/+y/-z/+z

[all_panels_ids, num_nodes] = lse_generate_panelIDs(idxS, L, M, N); 

%
% finding locations of surface panels (coordinates)
% This is needed for debug visualization (in the function) and for assigning the port nodes (later on)
% Note that 'all_panels_ids' is needed only for the debug visualization part
% 

all_panel_locs = lse_find_node_locations(num_nonair_cube, grid_intcon, L, M, N, dx, all_panels_ids);
  

%
% Finding the nodes on which excitation and ground points are defined
%

[nodeid_lft, nodeid_rght, nodeid_wlcond] = lse_find_port_nodes_from_positions(all_panel_locs, all_panels_ids, grid_intcon, L, M, N, dx, pnt_lft, pnt_rght, pnt_well_cond, fl_check_ports);  
 
 
% -------------------------------------------------------------------------
%                 Output file
% -------------------------------------------------------------------------

disp('------------------------------------------------------------------')
disp(['Generating output file ', fileoutname, ' ....']);

fout = fopen([fileoutname], 'w');

fprintf(fout, '* VoxHenry input file\n');
fprintf(fout, firstline);
fprintf(fout, '\n');

fprintf(fout, '* Frequency points (Hz)\n');
fprintf(fout, 'freq=');
for j = 1:size(freq_all,2)
    fprintf(fout, ' %g', freq_all(j));
    if ( mod(j,10) == 0 )
        fprintf(fout, '\n');
        fprintf(fout, 'freq=');
    end
end
fprintf(fout, '\n\n');

fprintf(fout, '* Voxel size (m)\n');
fprintf(fout, 'dx=%g\n', dx);
fprintf(fout, '\n');

fprintf(fout, '* Voxel grid dimension in voxel units; x, y, z\n');
fprintf(fout, 'LMN=%d,%d,%d\n', L, M, N);
fprintf(fout, '\n');

fprintf(fout, '* Voxel list\n');
fprintf(fout, '* Format is:\n');
fprintf(fout, '* V <index_x> <index_y> <index_z> <conductivity S/m>\n');
fprintf(fout, '*\n');
fprintf(fout, 'StartVoxelList\n');
% transform the linear intexes 'idxS' corresponding to non-empty voxels
% to voxel indices along x, y and z
[voxL,voxM,voxN] = ind2sub(size(sigma_e), idxS);  
for j = 1:size(idxS,1)
   % output the conductor voxel indices and conductivity
   fprintf(fout, 'V %d %d %d %g\r\n', voxL(j), voxM(j), voxN(j), sigma_e(voxL(j), voxM(j), voxN(j)));
end
fprintf(fout, 'EndVoxelList\n');

nodepos = ['-x'; '+x'; '-y'; '+y'; '-z'; '+z'];

fprintf(fout, '\n');
fprintf(fout, '* Port nodes list\n');
fprintf(fout, '* Format is:\n');
fprintf(fout, '* N <portname> <excitation or ground (P/N)> <voxel_index_x> <voxel_index_y> <voxel_index_z> <node (+z,-z,+x,-x,+y,-y)> \n');

for i_port = 1:num_ports
    fprintf(fout, '*\n');
    for j_node = 1:size(nodeid_lft{i_port},1)
        [voxL,voxM,voxN] = ind2sub(size(sigma_e), idxS(nodeid_lft{i_port}(j_node,2)));      
        fprintf(fout, 'N port%d P %d %d %d %s\n', i_port, voxL, voxM, voxN, nodepos(nodeid_lft{i_port}(j_node,3),:));
    end
    for j_node = 1:size(nodeid_rght{i_port},1)
        [voxL,voxM,voxN] = ind2sub(size(sigma_e), idxS(nodeid_rght{i_port}(j_node,2)));        
        fprintf(fout, 'N port%d N %d %d %d %s\n', i_port, voxL, voxM, voxN, nodepos(nodeid_rght{i_port}(j_node,3),:));
    end
end

% if there are grounded nodes
if (isempty(nodeid_wlcond) == 0)
    fprintf(fout, '\n');
    fprintf(fout, '* Grounded nodes list\n');
    fprintf(fout, '* Format is:\n');
    fprintf(fout, '* N <portname> <ground (G)> <voxel_index_x> <voxel_index_y> <voxel_index_z> <node (+z,-z,+x,-x,+y,-y)> \n');
    fprintf(fout, '*\n');

    for j_node = 1:size(nodeid_wlcond,1)
        [voxL,voxM,voxN] = ind2sub(size(sigma_e), idxS(nodeid_wlcond(:,2)));      
        fprintf(fout, 'N portGND G %d %d %d %s\n', voxL, voxM, voxN, nodepos(nodeid_wlcond(:,3)));
    end
end
 
fclose(fout);

disp(['Done!']);

