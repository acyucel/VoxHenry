function  K = volume_volume_sym(W,X,Y,Z,Xp,Yp,Zp,k0,r_m,r_n,dx,volume_ker,dir)

%% calculates the volume-volume integrals for the 3 integral types along x, y or z,
%% and returns a 3x1 matrix.
%% 'dir': 1 is x, 2 is y, 3 is z
%% Remark: in case of 'dir' along y or z, no need to recalculate pulse-pulse'
%%         therefore K(1) will be dummy

%% change of variables

% Jacobean
J = (dx/2)^6;

% sampling points for the Green kernel
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
        % apply l'hopital at the points where the kernel's evaluation is
        % 0/0
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


%% 3 unique kernels

K = zeros(3,1);

if(dir == 1)
    % change of variable for testing function
    Xt  = dx/2*X;

    % change of variable for basis function
    Xpb  = dx/2*Xp;

    % pulse-pulse'
    TB = 1;
    K(1) = sum(W .* TB .* ker) * J;

    % linear(x)-pulse'
    TB = Xt;
    K(2) = sum(W .* TB .* ker) * J;

    % linear(x)-linear(x')
    TB = Xpb.*Xt;
    K(3) = sum(W .* TB .* ker) * J;

elseif(dir == 2)
    % change of variable for testing function
    Yt  = dx/2*Y;

    % change of variable for basis function
    Ypb  = dx/2*Yp;
    
    % pulse-pulse' is dummy (assumed already calculated)
    K(1) = 0;

    % linear(y)-pulse'
    TB = Yt;
    K(2) = sum(W .* TB .* ker) * J;

    % linear(y)-linear(y')
    TB = Ypb.*Yt;
    K(3) = sum(W .* TB .* ker) * J;
    
else
    % change of variable for testing function
    Zt  = dx/2*Z;

    % change of variable for basis function
    Zpb  = dx/2*Zp;
    
    % pulse-pulse' is dummy (assumed already calculated)
    K(1) = 0;

    % linear(z)-pulse'
    TB = Zt;
    K(2) = sum(W .* TB .* ker) * J;

    % linear(z)-linear(z')
    TB = Zpb.*Zt;
    K(3) = sum(W .* TB .* ker) * J;
    
end

% end function
end
