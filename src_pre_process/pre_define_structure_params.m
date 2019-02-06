% Domain and unknown parameters
[L,M,N,~] = size(sigma_e); % domain size
nD = L*M*N; % number of variables in the system
% the voxels with non-zero conductivity
sigma_e_nonzero = abs(sigma_e(:)) > 1e-12;
if isempty(lambdaL)
    % the indices of the non-air voxel positions (where either sigma or lambdaL or both are non-null)
    idxS = find(sigma_e_nonzero);
else
    %the voxels with non-zero Londpn penetration depth 
    lambdaL_nonzero = abs(lambdaL(:)) > 1e-12; 
    % the indices of the non-air voxel positions (where either sigma or lambdaL or both are non-null)
    idxS = find(sigma_e_nonzero | lambdaL_nonzero); 
end
%idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions for 3 Cartesian components
idxS5 = [idxS; nD+idxS; 2*nD+idxS; 3*nD+idxS; 4*nD+idxS]; % for currents

% Constitutive parameters
% due to lowest frequency if multiple frequency is defined

% 'epsilon_r' 3D matrix is not strictly needed for the simulation - might use simply zero
% however it is used for post-processing visualization (through 'Mc', see VoxHenry_executer.m)
epsilon_r = ones(size(sigma_e));
epsilon_r(idxS) = 0.0;


