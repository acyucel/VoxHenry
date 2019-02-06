function [nodeid_4_grnd, nodeid_4_injectcurr]=lse_assign_exc_grnd_nodes(idxS, sigma_e, all_panels_ids, port_no, pnt_lft, pnt_rght, pnt_well_cond)

nodeid_4_injectcurr = [];
nodeid_4_grnd = [];

for jj=1:size(pnt_lft{port_no},1)
    % the 'port_no'-th cells of 'ptn_lft' contains the array of nodes belonging
    % to the 'port_no'-th port.
    % We scan each element of this array, that contains in the first three
    % columns the voxel position in L,M,N coordinates and in the fourth
    % column the node position w.r.t. the voxel, where 1 = -x, 2 = +x, 3 = -y etc.
    % and we convert the voxel postion to the non-empty voxel index of 'idxS';
    % based on that we find the actual panelID from 'all_panels_ids'
    voxelID = sub2ind(size(sigma_e), pnt_lft{port_no}(jj, 1), pnt_lft{port_no}(jj, 2), pnt_lft{port_no}(jj, 3));
    if(exist ('OCTAVE_VERSION', 'builtin') > 0)
        % as idxS are ordered, can lookup the voxelID in the non-empty voxel array 'idxS' quite fast (logN)
        nonemptyvoxelID = lookup(idxS, voxelID, 'm');
    else
        % this is linear search, as Matlab requires the 'Datafeed toolbox'
        % for using the function 'lookup' doing the binary search. So
        % resort to linear search.
        nonemptyvoxelID = find(idxS == voxelID);
        if isempty(nonemptyvoxelID)
            nonemptyvoxelID = 0;
        end
    end
    if nonemptyvoxelID == 0
        disp(['Internal error: cannot find left non-empty voxel corresponding to voxel ', num2str(pnt_lft{port_no}(jj, 1)), ',',  num2str(pnt_lft{port_no}(jj, 2)), ',', num2str(pnt_lft{port_no}(jj, 3))]);
    end
    nodeid = all_panels_ids(nonemptyvoxelID, pnt_lft{port_no}(jj,4));
    nodeid_4_injectcurr = [nodeid_4_injectcurr; nodeid];
end

for jj=1:size(pnt_rght{port_no},1)
    % the 'port_no'-th cells of 'ptn_rght' contains the array of nodes belonging
    % to the 'port_no'-th port.
    % We scan each element of this array, that contains in the first three
    % columns the voxel position in L,M,N coordinates and in the fourth
    % column the node position w.r.t. the voxel, where 1 = -x, 2 = +x, 3 = -y etc.
    % and we convert the voxel postion to the non-empty voxel index of 'idxS';
    % based on that we find the actual panelID from 'all_panels_ids'
    voxelID = sub2ind(size(sigma_e), pnt_rght{port_no}(jj, 1), pnt_rght{port_no}(jj, 2), pnt_rght{port_no}(jj, 3));
    if(exist ('OCTAVE_VERSION', 'builtin') > 0)
        % as idxS are ordered, can lookup the voxelID in the non-empty voxel array 'idxS' quite fast (logN)
        nonemptyvoxelID = lookup(idxS, voxelID, 'm');
    else
        % this is linear search, as Matlab requires the 'Datafeed toolbox'
        % for using the function 'lookup' doing the binary search. So
        % resort to linear search.
        nonemptyvoxelID = find(idxS == voxelID);
        if isempty(nonemptyvoxelID)
            nonemptyvoxelID = 0;
        end
    end
    if nonemptyvoxelID == 0
        disp(['Internal error: cannot find right non-empty voxel corresponding to voxel ', num2str(pnt_rght{port_no}(jj, 1)), ',',  num2str(pnt_rght{port_no}(jj, 2)), ',', num2str(pnt_rght{port_no}(jj, 3))]);
    end
    nodeid = all_panels_ids(nonemptyvoxelID, pnt_rght{port_no}(jj,4));
    nodeid_4_grnd = [nodeid_4_grnd; nodeid];
end

% this is instead a flat list of nodes to be grounded, not related to the ports
if (isempty(pnt_well_cond) == 0)
    % if not an empty list, however need to process only once, so we do it only
    % when 'port_no' is the first one (there must be at least one port)
    if port_no == 1
        for jj=1:size(pnt_well_cond,1)
            % We scan each element of the 'pnt_well_cond' array, that contains in the first three
            % columns the voxel position in L,M,N coordinates and in the fourth
            % column the node position w.r.t. the voxel, where 1 = -x, 2 = +x, 3 = -y etc.
            % and we convert the voxel postion to the non-empty voxel index of 'idxS';
            % based on that we find the actual panelID from 'all_panels_ids'
            voxelID = sub2ind(size(sigma_e), pnt_well_cond(jj, 1), pnt_well_cond(jj, 2), pnt_well_cond(jj, 3));
            if(exist ('OCTAVE_VERSION', 'builtin') > 0)
                % as idxS are ordered, can lookup the voxelID in the non-empty voxel array 'idxS' quite fast (logN)
                nonemptyvoxelID = lookup(idxS, voxelID, 'm');
            else
                % this is linear search, as Matlab requires the 'Datafeed toolbox'
                % for using the function 'lookup' doing the binary search. So
                % resort to linear search.
                nonemptyvoxelID = find(idxS == voxelID);
                if isempty(nonemptyvoxelID)
                    nonemptyvoxelID = 0;
                end
            end
            if nonemptyvoxelID == 0
                disp(['Internal error: cannot find non-empty voxel corresponding to voxel ', num2str(pnt_well_cond(jj, 1)), ',',  num2str(pnt_well_cond(jj, 2)), ',', num2str(pnt_well_cond(jj, 3))]);
            end
            nodeid = all_panels_ids(nonemptyvoxelID, pnt_well_cond(jj,4));
            % we accumulate these nodes to the ground nodes of the first port. No harm done by that
            % (all ground nodes of all ports are always kept to zero potential anyway)
            nodeid_4_grnd = [nodeid_4_grnd; nodeid]; 
        end
    end
end