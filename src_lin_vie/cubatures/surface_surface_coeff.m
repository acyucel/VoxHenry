function [I1_co, I2_co, I3_co, I4_co] = surface_surface_coeff(dx,k0)

%% multiplies the value of the non-constant result of the surface-surface integral with the respective constant coefficients

% normal vector components
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
            
for face = 1:6
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
       
       for l = 1:10
                    
                    % multiply by constant coefficient
                    I2_co(face,face_p,l) = coefficients(n,np,l,k0,dx,2);

                    % multiply by constant coefficient
                    I3_co(face,face_p,l) = coefficients(n,np,l,k0,dx,3);
                    
                    % multiply by constant coefficient
                    I4_co(face,face_p,l) = coefficients(n,np,l,k0,dx,4);

        end
        
        % multiply by constant coefficient
        I1_co(face,face_p) = coefficients(n,np,l,k0,dx,1);    

    end
    
end
  



end