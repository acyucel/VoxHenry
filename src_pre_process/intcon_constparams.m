function [idx,epsilon_r,sigma_e,grid_intcon] = intcon_constparams(r,Res,Cnt,Dims,Orients,e_r,s_e,fl_plot_intcons)
%%    Defines Constitutive parameters of interconnects
% _________________________________________________________________________
%
%       Defines interconnects in a given domain.
%       Optionally assigns constitutive parameters to interconnects
%
% _________________________________________________________________________
%
%% INPUT
%   r           4D (LxMxNx3) array with domain voxelized grid coordinates
%   Cnt         Cartesian coordinates of the centers of interconnects
%               row: interconnect id, column: dimensions
%   Dims        Dimensions of interconnects (LxWxH) (along x, y, and z)
%               row: interconnect id, column: length,width,height
%   Orients     Orientations of conductors, string
%   e_r         value of relative epsilon
%   s_e         value of electric conductivity
%
%% OUTPUT
%   idx         indices of the positions
%   epsilon_r   3D (LxMxN) array with relative epsilon
%   sigma_e     3D (LxMxN) array with electric conductivity
%
% -------------------------------------------------------------------------
%
%   Abdulkadir C. Yucel -- acyucel@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 4 )
    fprintf(1, '\n ERROR: not enough arguments\n');
    return
end
if(nargin < 5 || isempty(e_r))
    e_r = 1;
end
if(nargin < 6 || isempty(s_e))
    s_e = 0;
end

% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

[L,M,N,~] = size(r);

num_intcon=size(Cnt,1);

% define bounds for each interconnect
x_bnd=zeros(num_intcon,2);y_bnd=zeros(num_intcon,2);z_bnd=zeros(num_intcon,2);
for kk=1:num_intcon
    temp_cen=Cnt(kk,1:3);
    temp_dim=Dims(kk,1:3);
    
    if (Orients(kk) == 'x')
        x_bnd(kk,1:2)=[temp_cen(1)-temp_dim(1)*0.5 temp_cen(1)+temp_dim(1)*0.5];
        y_bnd(kk,1:2)=[temp_cen(2)-temp_dim(2)*0.5 temp_cen(2)+temp_dim(2)*0.5];
        z_bnd(kk,1:2)=[temp_cen(3)-temp_dim(3)*0.5 temp_cen(3)+temp_dim(3)*0.5];
    elseif (Orients(kk) == 'y')
        x_bnd(kk,1:2)=[temp_cen(1)-temp_dim(2)*0.5 temp_cen(1)+temp_dim(2)*0.5];
        y_bnd(kk,1:2)=[temp_cen(2)-temp_dim(1)*0.5 temp_cen(2)+temp_dim(1)*0.5];
        z_bnd(kk,1:2)=[temp_cen(3)-temp_dim(3)*0.5 temp_cen(3)+temp_dim(3)*0.5];
    elseif (Orients(kk) == 'z')
        x_bnd(kk,1:2)=[temp_cen(1)-temp_dim(3)*0.5 temp_cen(1)+temp_dim(3)*0.5];
        y_bnd(kk,1:2)=[temp_cen(2)-temp_dim(2)*0.5 temp_cen(2)+temp_dim(2)*0.5];
        z_bnd(kk,1:2)=[temp_cen(3)-temp_dim(1)*0.5 temp_cen(3)+temp_dim(1)*0.5];
    else
        error('Orients should have only x, y, or z strings')
    end
    
end

boolean_tens=zeros(L,M,N);

tola=1e-12;
for ll=1:L
    for mm=1:M
        for nn=1:N
            for kk=1:num_intcon
                
                if ( r(ll,mm,nn,1) > x_bnd(kk,1)-tola && r(ll,mm,nn,1) < x_bnd(kk,2)+tola &&...
                        r(ll,mm,nn,2) > y_bnd(kk,1)-tola && r(ll,mm,nn,2) < y_bnd(kk,2)+tola && ...
                        r(ll,mm,nn,3) > z_bnd(kk,1)-tola && r(ll,mm,nn,3) < z_bnd(kk,2)+tola)
                    
                    boolean_tens(ll,mm,nn)=1;
                    
                end
            end
        end
    end
end

idx = find(boolean_tens); % get indices of elements

% -------------------------------------------------------------------------
% Assign output
% -------------------------------------------------------------------------

epsilon_r = ones(L,M,N);
epsilon_r(idx) = e_r;

sigma_e = zeros(L,M,N);
sigma_e(idx) = s_e;

% -------------------------------------------------------------------------
% Assign the grid of conductors and visualize it 
% -------------------------------------------------------------------------

[L,M,N,~] = size(r);
grid_intcon = zeros(L,M,N,3);
grid_intcon(idx) = r(idx);
grid_intcon(L*M*N+idx) = r(L*M*N+idx);
grid_intcon(2*L*M*N+idx) = r(2*L*M*N+idx);

% disp('Attention ::: Pad NaNs to grid_intcon for air voxels instead zeros! ')
% disp('This will avoid possible bug that appears when you have a voxel centered at (0,0,0)')

% plot geometry
%fl_plot_intcons=1;
if (fl_plot_intcons == 1)
    tic
    plot_boxes_of_grid(grid_intcon,Res);
    disp(['Time for visualizing the structure ::: ',num2str(toc)])
end


