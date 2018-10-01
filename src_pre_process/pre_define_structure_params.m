
% Constitutive parameters
% due to lowest frequency if multiple frequency is defined
Mr = epsilon_r - 1j*sigma_e/(eo*omega); % permittivity
Mc = Mr - 1.0; % susceptibility
OneoverMc = 1.0 ./ Mc; % one over susceptibility

% Domain and unknown parameters
[L,M,N,~] = size(r); % domain size
nD = L*M*N; % number of variables in the system
dx = Res; % voxel size
idxS = find(abs(Mc(:)) > 1e-12); % the indices of the non-air voxel positions
%idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions for 3 Cartesian components
idxS5 = [idxS; nD+idxS; 2*nD+idxS; 3*nD+idxS; 4*nD+idxS]; % for currents