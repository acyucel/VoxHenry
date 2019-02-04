function [sigma_e, lambdaL, freq, dx, num_ports, pnt_lft, pnt_rght, pnt_well_cond] = pre_input_file(fileinname)

tini = tic;

% constants

% average line length
averlinelen = 22;
% sides
sidetypes = {'-x', '+x', '-y', '+y', '-z', '+z'};

% open file and get info
[fid,errmsg] = fopen(fileinname, 'r');
if fid < 0
    disp(['ERROR opening the input file ', fileinname]);
    disp(['Error is: ', errmsg]);
    return;
end
  
% init 
freq = [];
dx = 1e-6;
pnt_lft = {};
pnt_rght = {};
pnt_well_cond = [];
portnames = {};
supercond = 0;

% get first line 
line = fgetl(fid);
% testing ischar() instead of feof(), as fgetl() and fgets() sometimes set 
% the end-of-file indicator before they return a value of -1.
% (as described in MatLab docs)
% Note that this behavior does not conform to the ANSI specifications for
% the related C language functions.
while ischar(line) == 1
    % split the line on the tokens (after removing trailing and leading spaces)
    % note that we include '=' in the delimiters
    scanline = strsplit(strtrim(line),  {' ','\f','\n','\r','\t','\v','=',','});
    if ~isempty(scanline)  
        % if frequency definitions
        if strcmpi(scanline{1},'freq') == 1
            for index = 2:size(scanline,2)
                % append the newly read frequency
                freq = [freq str2double(scanline{index})];
            end
  	 	elseif strcmpi(scanline{1},'dx') == 1
            if size(scanline,2) >= 2
                dx = str2double(scanline{2});
            end
  	 	elseif strcmpi(scanline{1},'LMN') == 1
          if size(scanline,2) >= 4
                Lsize = str2double(scanline{2});
                Msize = str2double(scanline{3});
                Nsize = str2double(scanline{4});
                % pre-allocate sigma_e
                sigma_e = zeros(Lsize, Msize, Nsize);
                % pre-allocate lambdaL
                lambdaL = zeros(Lsize, Msize, Nsize);
          end
  	 	elseif strcmpi(scanline{1},'superconductor') == 1
          supercond = 1;
  	 	elseif strcmpi(scanline{1},'startvoxellist') == 1
            % fast scan of voxel list
            if supercond == 0
                % flags mean: do not stop on error
                voxellist = textscan(fid, 'V %d %d %d %f', 'ReturnOnError', 1);
                sigma_e(sub2ind(size(sigma_e), voxellist{1}, voxellist{2}, voxellist{3})) = voxellist{4};
            else
                % flags mean: do not stop on error
                voxellist = textscan(fid, 'V %d %d %d %f %f', 'ReturnOnError', 1);
                sigma_e(sub2ind(size(sigma_e), voxellist{1}, voxellist{2}, voxellist{3})) = voxellist{4};
                % lambda penetration depth in case of superconductors
                lambdaL(sub2ind(size(lambdaL), voxellist{1}, voxellist{2}, voxellist{3})) = voxellist{5};
            end
  	 	elseif strcmpi(scanline{1},'N') == 1
            if size(scanline,2) >= 7
                pname = scanline{2};
                % check the port number based on its name
                portnumber = find(strcmpi(portnames, pname));
                % if this is a new port name
                if isempty(portnumber)
                    % add it to the list of port names
                    portnames = {portnames{:} pname};
                    portnumber = size(portnames,2);
                    % create corresponding new cells in the left and right cell arrays
                    pnt_lft{portnumber} = [];
                    pnt_rght{portnumber} = [];
                end
                xindex = sscanf(scanline{4}, '%d');
                yindex = sscanf(scanline{5}, '%d');
                zindex = sscanf(scanline{6}, '%d');
                side = find(strcmp(sidetypes, scanline{7}));
                if strcmpi(scanline{3}, 'P') == 1
                    % add the node to the list of nodes belonging to this port
                    pnt_lft{portnumber} = [pnt_lft{portnumber}; xindex yindex zindex side];
                elseif strcmpi(scanline{3}, 'N') == 1
                    % add the node to the list of nodes belonging to this port
                    pnt_rght{portnumber} = [pnt_rght{portnumber}; xindex yindex zindex side];
                else
                    % add the node to the list of nodes to be grounded
                    pnt_well_cond = [pnt_well_cond; xindex yindex zindex side];
                end
            end
        end
    end
   
    line = fgetl(fid);
end

% if no superconductor, let's empty the 'lambdaL' tensor, so we know these are standard conductors 
if supercond == 0
    lambdaL = []
end

% sanity checks
if size(pnt_lft,2) ~= size(pnt_rght,2)
    disp(['Input file error: number of left ports is ', num2str(size(pnt_lft,2)), ' while number of right ports is ', num2str(size(pnt_rght,2))]);
end
num_ports = size(pnt_lft,2);

if isempty(freq)
    disp('Input file error: no frequencies specified. Defaulting to 1 Hz');
    freg = 1;
end

fclose(fid);

tend = toc(tini);
disp(['Time for reading input file::: ' ,num2str(tend)]);
disp('-----------------------------------------------------')


