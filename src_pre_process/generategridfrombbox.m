function [r] = generategridfrombbox(res,x,y,z,fl_plot_boxes_w_grid)
%    Generates a 3D grid from the limits of a given domain
% _________________________________________________________________________
%
%       Generates a 3D cartesian grid of coordinates with given resolution
%       The resolution is fixed, and the final domain is the minimum for
%       the given resolution that encloses the specified dimensions
%
% _________________________________________________________________________
%
%% INPUT
%   res         resolution
%   x           minimum and maximum values of x
%   y           minimum and maximum values of y
%   z           minimum and maximum values of z
%
%
%% OUTPUT
%   r           4D (LxMxNx3) array with domain voxelized grid coordinates
%
%
% -------------------------------------------------------------------------
%
%   Abdulkadir C. Yucel -- acyucel@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________

%clc;close all;clear all
%res = 1e-3;
%x=[0 4e-3];y=[-0.5e-3 0.5e-3];z=[-0.5e-3 0.5e-3];

% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 2 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 3 || isempty(y))
   y = x;
end
if(nargin < 4 || isempty(z))
   z = x;
end
if(nargin < 5 || isempty(fl_plot_boxes_w_grid))
   fl_plot_boxes_w_grid = 1;
end

% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

% just in case
x = squeeze(x);
y = squeeze(y);
z = squeeze(z);

% if the input is a single element, takes negative and positive
if length(x) == 1
    x = [-x x];
end
if length(y) == 1
    y = [-y y];
end
if length(z) == 1
    z = [-z z];
end

% -------------------------------------------------------------------------
% Obtain bounding box
% -------------------------------------------------------------------------

% 1) Computing edge lengths of bounding box
bbox_elem(1)=x(2)-x(1);
bbox_elem(2)=y(2)-y(1);
bbox_elem(3)=z(2)-z(1);

%2) finding number of boxes along x, y, and z axes

for kk=1:3
    nelem_xyz(kk)=floor(bbox_elem(kk)/res);
end
   
% 3) Print out the location of local origin
bbox_origin=[min(x) min(y) min(z)];
disp(['Location of local origin ::  ', num2str(bbox_origin,'%10.5e %10.5e %10.5e')])
disp(['# of boxes along x, y, and z directions::  ', num2str(nelem_xyz,'%10i %10i %10i')])

% 4) Finding the centers of boxes along all directions

for kk=1:nelem_xyz(1)
    xx(kk)=(((kk-1)*res+kk*res)*0.5)+bbox_origin(1);
end
for kk=1:nelem_xyz(2)
    yy(kk)=(((kk-1)*res+kk*res)*0.5)+bbox_origin(2);
end
for kk=1:nelem_xyz(3)
    zz(kk)=(((kk-1)*res+kk*res)*0.5)+bbox_origin(3);
end

% -------------------------------------------------------------------------
% Generate grid
% -------------------------------------------------------------------------

r = grid3d(xx,yy,zz);


% -------------------------------------------------------------------------
% Plot grid with boxes
% -------------------------------------------------------------------------

%fl_plot_boxes_w_grid=0;

if (fl_plot_boxes_w_grid == 1)
    tic
    plot_boxes_of_grid(r,res);
    disp(['Time for visualizing the computational domain ::: ',num2str(toc)])
end

