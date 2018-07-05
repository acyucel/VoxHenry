function [W,X,Y,Z,Xp,Yp,Zp] = weights_points(N,dim)
%%
% dim: dimension of integration
% N:   number of points for one dimension
% returns the N^dim points according to Gauss-Legendre quadrature (X,Y,Z,X',Y',Z')
% and the associated weights W


%% Clenshaw - Curtis
% % weights
% w = zeros(N,1);
% w(1) = 2;
% w(3:2:N) = 2./(1-(2:2:N-1).^2);
% q = [w;w(N-1:-1:2)];
% w1 = real(ifft(q));
% w1d = [w1(1); 2*w1(2:N-1); w1(N)];
% 
% % points
% t = pi/(N-1)*(0:N-1)';
% x1d = cos(t);


%% Gauss-Legendre
[ w1d , x1d ] = gauss_1d(N);
% points
x1d = transpose(x1d);
% weights
w1d = transpose(w1d);


%% expansion up to 6 dimensions
if dim == 1
    
    W = w1d;
    X = x1d;
    
elseif dim == 2
    
    w2d = w1d*transpose(w1d);
    W = w2d(:);
    [X,Y] = meshgrid(x1d);
    X = X(:);
    Y = Y(:);
    
elseif dim == 3
    
    w2d = w1d*transpose(w1d);
    w3d = w2d(:)*transpose(w1d);
    W = w3d(:);
    [X,Y,Z] = meshgrid(x1d);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    
elseif dim == 4
    
    w2d = w1d*transpose(w1d);
    w3d = w2d(:)*transpose(w1d);
    w4d = w3d(:)*transpose(w1d);
    W = w4d(:);
    [X,Y,Z,Xp] = ndgrid(x1d);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    Xp = Xp(:);  
    
elseif dim == 5
    
    w2d = w1d*transpose(w1d);
    w3d = w2d(:)*transpose(w1d);
    w4d = w3d(:)*transpose(w1d);
    w5d = w4d(:)*transpose(w1d);
    W = w5d(:);
    [X,Y,Z,Xp,Yp] = ndgrid(x1d);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    Xp = Xp(:);
    Yp = Yp(:);
    
elseif dim == 6
    
    w2d = w1d*transpose(w1d);
    w3d = w2d(:)*transpose(w1d);
    w4d = w3d(:)*transpose(w1d);
    w5d = w4d(:)*transpose(w1d);
    w6d = w5d(:)*transpose(w1d);
    W = w6d(:);
    [X,Y,Z,Xp,Yp,Zp] = ndgrid(x1d);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    Xp = Xp(:);
    Yp = Yp(:);
    Zp = Zp(:);
    
end
    
    
end