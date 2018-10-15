function I_S = surface_surface_kernels(I1_co,I2_co,I3_co,I4_co,W,A,B,C,D,Np,k0,r_m,r_n,m,dx,volume_ker)

% calculates the near interactions for the 10 integral types by decomposing
% the volume-volume integrals to surface-surface ones and accordingly
% calculate the singular ones with directfn by calling mexed functions

% I1_co,I2_co,I3_co,I4_co: constant coefficients
% volumer_ker: the corresponding kernels of table 1,2,3. 1: G, 2: 1/R, 3: R


% 4D integral jacobian
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

% scan the observation voxel faces
for face = 1:6
    % scane the source voxel faces
    for face_p = 1:6
        
        % get the oredered vertices of the squares for the specific
        % square-square adjacency type
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

          
        
        
        % scan over the integral kernels
        for l = 1:10
            
            %% 1-st Kernel (1st row of the tables 1,2,3)
            
            % if the coefficient infront of the integral is nonzero
            % caculate the surface-surface kernel
            if I1_co(face,face_p) ~= 0
                
                % if the integrals is non-singular use regular integration routines.
                % Not all the 36 face-face interaction are singular. They
                % are singular only when the 2 faces share a common
                % vertice, edge or face.
                if adjacency_type(face,face_p) == 0

                    % evaluate kernel at the integration points
                    f1 = kernels(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,n,l,1,volume_ker);
                    I1(face,face_p,l) = I1_co(face,face_p) * sum( W .*f1 ) * J;
                
                % in essence the singular surface-surface integrals are the
                % slow calculations...
                    
                % self term singular integral (common face)
                elseif adjacency_type(face,face_p) == 1
                    % NOTE the l-1 indexing which is needed because of the
                    % indexing difference between c++ and matlab.
                    I1(face,face_p,l) = I1_co(face,face_p) * singular_ST(volume_ker,l-1,k0,dx,Np,1,n,np,R(:,1),R(:,2),R(:,3),R(:,4));

                % edge adjacent singular integral (common edge)
                elseif adjacency_type(face,face_p) == 2
                    
                    I1(face,face_p,l) = I1_co(face,face_p) * singular_EA(volume_ker,l-1,k0,dx,Np,1,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6));

                % vertex adjacent singular integral ( common vertex)   
                elseif adjacency_type(face,face_p) == 3
                    I1(face,face_p,l) = I1_co(face,face_p) * singular_VA(volume_ker,l-1,k0,dx,Np,1,n,np,R(:,1),R(:,2),R(:,3),R(:,4),R(:,5),R(:,6),R(:,7));

                end
                
            end
            
            
            %% 2-nd Kernel (2nd row of the tables 1,2,3)
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
            
            
            %% 3-rd Kernel (3rd row of the tables 1,2,3)
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
            
            
            
        %% 4-th Kernel (4th row of the tables 1,2,3)
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
                
            
        end
    end
end

% take the sum of the 36 face-face integrals according to reports eq. 15 
I_SK = zeros(10,1);
for l = 1:10
    
    I_SK(l) = sum(sum(I1(:,:,l))) + sum(sum(I2(:,:,l))) + sum(sum(I3(:,:,l))) + sum(sum(I4(:,:,l)));
    
end

% calculate the integrals
I_S(1) =  I_SK(1);
I_S(2) =  I_SK(2);
I_S(3) = -I_SK(4);
I_S(4) =  I_SK(8) + I_SK(9);
I_S(5) = -2*I_SK(6);
I_S(6) =  I_SK(8) - I_SK(9);
I_S(7)=  I_SK(8) + I_SK(9) + 4*I_SK(10);

end
