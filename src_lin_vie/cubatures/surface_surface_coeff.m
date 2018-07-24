function [I1_co, I2_co, I3_co, I4_co] = surface_surface_coeff(dx,k0)

% I1_co = -n * n'
% I2_co = n' * \nabla f_m 
% I3_co = n' * \nabla' f_n 
% I4_co = (n' * \nabla f_m) * ( n' * \nabla' f_n )

%% calculates the constant coefficients infront of 4 surface-surface integral kernels

% normal vector components of each face of the cube
nx = [-1 1  0 0  0 0];
ny = [ 0 0 -1 1  0 0];
nz = [ 0 0  0 0 -1 1];


% n,n'
I1_co = zeros(6,6);
% n,n',type
I2_co = zeros(6,6,10);
% n,n',type
I3_co = zeros(6,6,10);
% n,n',type
I4_co = zeros(6,6,10);

% scan the faces of the observation voxel
for face = 1:6
    % scan the faces of the source voxel
    for face_p = 1:6
        
       n = zeros(3,1);
       np = zeros(3,1);
        
       % calculate the normal vectors 
       n (1)  = nx(face);
       np(1)  = nx(face_p);
       n (2)  = ny(face);
       np(2)  = ny(face_p);
       n (3)  = nz(face);
       np(3)  = nz(face_p);
       
       % scane the products between the basis and testing functions
       for l = 1:10
                    
                    % calculate constant coefficients for second
                    % surface-surface kernel
                    I2_co(face,face_p,l) = coefficients(n,np,l,k0,dx,2);

                    % calculate constant coefficients for third kernel
                    I3_co(face,face_p,l) = coefficients(n,np,l,k0,dx,3);
                    
                    % calculate constant coefficients for fourth kernel
                    I4_co(face,face_p,l) = coefficients(n,np,l,k0,dx,4);

        end
        
        % calculate constant coefficients for first kernel
        I1_co(face,face_p) = coefficients(n,np,l,k0,dx,1);    

    end
    
end
  



end