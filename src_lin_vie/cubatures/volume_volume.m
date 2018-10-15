function  I = volume_volume(W,X,Y,Z,Xp,Yp,Zp,k0,r_m,r_n,dx,volume_ker)

%% calculates the volume-volume integrals for the 6  integral types and returns a 6x1 matrix

%% change of variables

% Jacobian
J = (dx/2)^6;

% sampling points for the Green kernel after change of variable
Xg  = r_m(1) + dx/2*X;
Yg  = r_m(2) + dx/2*Y;
Zg  = r_m(3) + dx/2*Z;

Xpg = r_n(1) + dx/2*Xp;
Ypg = r_n(2) + dx/2*Yp;
Zpg = r_n(3) + dx/2*Zp;

% sampling the functions at the appropriate points and change of varibale to transform the integration domain to unit cube-cube
R = sqrt( (Xg-Xpg).^2 + (Yg-Ypg).^2 + (Zg-Zpg).^2 );

switch volume_ker
    case 1
        % kernel for far integration
        ker = 1/4/pi * exp(-1j*k0*R)./ R;
    case 2
        % 1st subtracted kernel
        ker = 1/4/pi ./ R;
    case 3
        % 2nd subtracted kernel
        ker = (-k0^2/8/pi) * R;
    case 4
        % smoothed kernel
        ker1 = 1/4/pi * ( exp(-1j*k0*R) - 1.0 ) ./R;
        % apply l'hopital at the points where the kernel's evaluation is 0/0
        ker1(isnan(abs(ker1))) = -1j*k0 /4/pi;
        
        ker2 = k0^2/8/pi * R;
        ker = ker1 + ker2;
end


% change of variable for testing function
Xt  = dx/2*X;
Yt  = dx/2*Y;
Zt  = dx/2*Z;

% change of variable for basis function
Xpb  = dx/2*Xp;
Ypb  = dx/2*Yp;
Zpb  = dx/2*Zp;


%% 10 unique kernels

% pulse-pulse'
TB = 1;
K(1) = sum(W .* TB .* ker) * J;

% Exploiting integral symmetry. K(2) =  -K(3)
% pulse-linear(x')
%TB = Xpb;
%K(2) = sum(W .* TB .* ker) * J;

% linear(x)-pulse'
TB = Xt;
K(3) = sum(W .* TB .* ker) * J;

% Exploiting integral symmetry. K(4) =  -K(5)
% pulse-linear(y')
%TB = Ypb;
%K(4) = sum(W .* TB .* ker) * J;

% linear(y)-pulse'
TB = Yt;
K(5) = sum(W .* TB .* ker) * J;

% Exploiting integral symmetry. K(6) =  -K(7)
% pulse-linear(z')
%TB = Zpb;
%K(6) = sum(W .* TB .* ker) * J;

% linear(z)-pulse'
TB = Zt;
K(7) = sum(W .* TB .* ker) * J;

% linear(x)-linear(x')
TB = Xpb.*Xt;
K(8) = sum(W .* TB .* ker) * J;

% linear(y)-linear(y')
TB = Ypb.*Yt;
K(9) = sum(W .* TB .* ker) * J;

% linear(z)-linear(z')
TB = Zpb.*Zt;
K(10) = sum(W .* TB .* ker) * J;


% calculate the system's integrals
%
% Gx,x
I(1) =  K(1);
% G2D,x
I(2) =  -K(3);
% G2D,y
I(3) = K(5);
% G2D,2D
I(4) =  K(8) + K(9);
% G3D,z
I(5) = 2*K(7);
% G2D,3D
I(6) =  K(8) - K(9);
% G3D,3D
I(7)=  K(8) + K(9) + 4*K(10);

end
