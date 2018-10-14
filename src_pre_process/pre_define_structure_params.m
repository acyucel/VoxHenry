% Domain and unknown parameters
[L,M,N,~] = size(sigma_e); % domain size
nD = L*M*N; % number of variables in the system
idxS = find(abs(sigma_e(:)) > 1e-12); % the indices of the non-air voxel positions
%idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions for 3 Cartesian components
idxS5 = [idxS; nD+idxS; 2*nD+idxS; 3*nD+idxS; 4*nD+idxS]; % for currents

% Constitutive parameters
% due to lowest frequency if multiple frequency is defined

% 'epsilon_r' 3D matrix is not strictly needed - might use simply zero
epsilon_r = ones(size(sigma_e));
epsilon_r(idxS) = 0.0;


