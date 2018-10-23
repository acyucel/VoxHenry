function K = surface_surface_kernels_sym(I1_co,I2_co,I3_co,I4_co,W,A,B,C,D,Np,k0,r_m,r_n,m,dx,volume_ker,dir)

%% calculates the surface-surface integrals for the 6 integral types for the
%  4 surface-surface kernels and returns 3 6x6x10 matrices and 1 6x6 for the
%  36 faces interactions where the last integral is separated in its 2
%  different kernels. The constant coefficients in front of each kernel
%  are not taken into account here, this is doen later on
%% 'dir': 1 is x, 2 is y, 3 is z

% 3 unique kernels
%
if(dir == 1)
    % along x, we are just interested in:
    % - pulse-pulse'          ->   I_SK(1)
    % - linear(x)-pulse'      ->   I_SK(3)
    % - linear(x)-linear(x')  ->   I_SK(8)
    l_vals = [1 3 8];
 elseif(dir ==2)
    % along y, we are just interested in:
    % pulse-pulse' is dummy (assumed already calculated)
    % - linear(y)-pulse'      ->   I_SK(5)
    % - linear(y)-linear(y')  ->   I_SK(9)
    K(1) = 0;
    l_vals = [5 9];
else 
    % along y, we are just interested in:
    % - linear(z)-pulse'      ->   I_SK(7)
    % - linear(z)-linear(z')  ->   I_SK(10)
    K(1) = 0;
    l_vals = [7 10];
end

% 4D jacobean
J = (dx/2)^4;

% normal vector components
nx = [-1 1  0 0  0 0];
ny = [ 0 0 -1 1  0 0];
nz = [ 0 0  0 0 -1 1];

% n,n',type
I1 = zeros(6,6,10);
% n,n',type
I2 = zeros(6,6,10);
% n,n',type
I3 = zeros(6,6,10);
% n,n',type
I4 = zeros(6,6,10);


% calculate the adjacency type (ST,EA,VA) for all the square-square
% interactions of the faces of the voxel-voxel and the (4,6,7) points in  
% proper order (as required at directfn) for the calculation of the singular integrals
[ adjacency_type , ordered_points ] = points_mappping(m,dx);

for face = 1:6
    for face_p = 1:6
        
        % get the points for the specific square-square
        R = ordered_points(:,:,face,face_p);
        
        n  = zeros(3,1);
        np = zeros(3,1);
        
        % calculate the normal vectors
        n (1)  = nx(face);
        np(1)  = nx(face_p);
        n (2)  = ny(face);
        np(2)  = ny(face_p);
        n (3)  = nz(face);
        np(3)  = nz(face_p);
        
        % calculate the centre of the squares
        r_ms = r_m + n *dx/2;
        r_ns = r_n + np*dx/2;
        
        % calculate the points for all 6 variables
        [X,Y,Z,Xp,Yp,Zp] = points_const_4D(A,B,C,D,r_ms,r_ns,dx,n,np); 
        
        for l = l_vals
            
            %% 1-st Kernel
            if I1_co(face,face_p) ~= 0
                
                % non-singular square-square 
                if adjacency_type(face,face_p) == 0

                    % evaluate kernel at the points of interest
                    f1 = kernels(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,n,l,1,volume_ker);
                    I1(face,face_p,l) = I1_co(face,face_p) * sum( W .*f1 ) * J;

                % singular self term
                elseif adjacency_type(face,face_p) == 1
                   
                    I1(face,face_p,l) = I1_co(face,face_p) * singular_ST(volume_ker,l-1,k0,dx,Np,1,n,np,R(:,1),R(:,2),R(:,3),R(:,4));

                % singular edge adjacent
                elseif adjacency_type(face,face_p) == 2
                    
                    I1(face,face_p,l) = I1_co(face,face_p) * singular_EA(volume_ker,l-1,k0,dx,Np,1,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6));

                % singular vertex adjacent    
                elseif adjacency_type(face,face_p) == 3
                    
                    I1(face,face_p,l) = I1_co(face,face_p) * singular_VA(volume_ker,l-1,k0,dx,Np,1,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7));

                end
                
            end
            
            
            %% 2-nd Kernel
            if I2_co(face,face_p,l) ~= 0
                
                % non-singular square-square 
                if adjacency_type(face,face_p) == 0

                    % evaluate kernel at the points of interest
                    f2 = kernels(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,n,l,2,volume_ker);
                    I2(face,face_p,l) = I2_co(face,face_p,l) * sum( W .*f2 ) * J;

                % singular self term
                elseif adjacency_type(face,face_p) == 1
                   
                    I2(face,face_p,l) = I2_co(face,face_p,l) * singular_ST(volume_ker,l-1,k0,dx,Np,2,n,np,R(:,1),R(:,2),R(:,3),R(:,4));

                % singular edge adjacent
                elseif adjacency_type(face,face_p) == 2
                    
                    I2(face,face_p,l) = I2_co(face,face_p,l) * singular_EA(volume_ker,l-1,k0,dx,Np,2,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6));

                % singular vertex adjacent    
                elseif adjacency_type(face,face_p) == 3
                    
                    I2(face,face_p,l) = I2_co(face,face_p,l) * singular_VA(volume_ker,l-1,k0,dx,Np,2,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7));

                end
                
            end
            
            
            %% 3-rd Kernel
            if I3_co(face,face_p,l) ~= 0
                
                % non-singular square-square 
                if adjacency_type(face,face_p) == 0

                    % evaluate kernel at the points of interest
                    f3 = kernels(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,n,l,3,volume_ker);
                    I3(face,face_p,l) = I3_co(face,face_p,l) * sum( W .*f3 ) * J;

                % singular self term
                elseif adjacency_type(face,face_p) == 1
                   
                    I3(face,face_p,l) = I3_co(face,face_p,l) * singular_ST(volume_ker,l-1,k0,dx,Np,3,n,np,R(:,1),R(:,2),R(:,3),R(:,4));

                % singular edge adjacent
                elseif adjacency_type(face,face_p) == 2
                    
                    I3(face,face_p,l) = I3_co(face,face_p,l) * singular_EA(volume_ker,l-1,k0,dx,Np,3,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6));

                % singular vertex adjacent    
                elseif adjacency_type(face,face_p) == 3
                    
                    I3(face,face_p,l) = I3_co(face,face_p,l) * singular_VA(volume_ker,l-1,k0,dx,Np,3,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7));

                end
                
            end
            
            %% 4-th Kernel 
            if I4_co(face,face_p,l) ~= 0 
                
                % non-singular square-square         
                if adjacency_type(face,face_p) == 0

                    % evaluate kernel at the points of interest
                    f4 = kernels(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,n,0,4,volume_ker);
                    I4(face,face_p,l) = I4_co(face,face_p,l) * sum( W .*f4 ) * J;
                    
                % singular self term
                elseif adjacency_type(face,face_p) == 1
                    
                    I4(face,face_p,l)   = I4_co(face,face_p,l) * singular_ST(volume_ker,0,k0,dx,Np,4,n,np,R(:,1),R(:,2),R(:,3),R(:,4));
                    
                % singular edge adjacent
                elseif adjacency_type(face,face_p) == 2
                    
                    I4(face,face_p,l)   = I4_co(face,face_p,l) * singular_EA(volume_ker,0,k0,dx,Np,4,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6));
                    
                % singular vertex adjacent    
                elseif adjacency_type(face,face_p) == 3
                    
                    I4(face,face_p,l)   = I4_co(face,face_p,l) * singular_VA(volume_ker,0,k0,dx,Np,4,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7));
                
                end
            
            end     

        end %for
    end
end

% take the sum of the 36 face-face integrals according to reports eq. 15 
I_SK = zeros(10,1);
for l = l_vals
    
    I_SK(l) = sum(sum(I1(:,:,l))) + sum(sum(I2(:,:,l))) + sum(sum(I3(:,:,l))) + sum(sum(I4(:,:,l)));
    
end

if(dir == 1)
  K = I_SK(l_vals);
else
  K = zeros(3,1);
  K(2:3) = I_SK(l_vals);
end


end
