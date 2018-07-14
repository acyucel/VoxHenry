function f = kernels(X,Y,Z,Xp,Yp,Zp,k0,dx,r_m,r_n,n,l,ker_type,volume_ker)

%% calculates the surface-surface integrals for the 4 kernels for each testing and basis functions combination 

% chose one of the four surface-surface kernels (rows of the tables 1,2,3)
switch ker_type
    case 1
        %% 1st surface-surface kernel
        switch l
            case 1
                fm = 1;
                fn = 1;
            case 2
                fm = 1;
                fn = Xp - r_n(1);
            case 3
                fm = X  - r_m(1); 
                fn = 1;
            case 4
                fm = 1;
                fn = Yp - r_n(2);
            case 5
                fm = Y  - r_m(2);
                fn = 1;
            case 6
                fm = 1;
                fn = Zp - r_n(3);
            case 7 
                fm = Z  - r_m(3);
                fn = 1;
            case 8
                fm = X  - r_m(1); 
                fn = Xp - r_n(1);
            case 9
                fm = Y  - r_m(2);
                fn = Yp - r_n(2);
            case 10
                fm = Z  - r_m(3);
                fn = Zp - r_n(3);
        end
        
        switch volume_ker
            % reduced kernels from G
            case 1
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                G  = 1/4/pi * exp(-1i*k0*R) ./ R;
                G0 = 1/4/pi ./ R;
                f = fm.*fn .* (G-G0) / (1j*k0)^2;
            % reduced kernels from 1/R  
            case 2
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                f = fm.*fn .* R/2 /4/pi;
            % reduced kernels from R   
            case 3
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                f = fm.*fn .* R.^3/12 *(-k0^2/8/pi);
        end
        
    case 2
        %% 2nd surface-surface kernel
        % fn
        switch l
            case 1
                fn = 1;
            case 2
                fn = Xp - r_n(1);
            case 3
                fn = 1;
            case 4
                fn = Yp - r_n(2);
            case 5
                fn = 1;
            case 6
                fn = Zp - r_n(3);
            case 7 
                fn = 1;
            case 8
                fn = Xp - r_n(1);
            case 9
                fn = Yp - r_n(2);
            case 10
                fn = Zp - r_n(3);
        end
            
        
        switch volume_ker
            % reduced kernels from G
            case 1
                % Green
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                G = 1/4/pi * exp(-1i*k0*R) ./ R;
                G0 = 1/4/pi ./ R;

                % F
                ff(:,1) = (X-Xp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
                ff(:,2) = (Y-Yp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
                ff(:,3) = (Z-Zp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);

                % F - RG0/2
                gg(:,1) = (X-Xp)/2.*G0;
                gg(:,2) = (Y-Yp)/2.*G0;
                gg(:,3) = (Z-Zp)/2.*G0;
                hh = ff - gg;
                q = n(1)*hh(:,1) + n(2)*hh(:,2) + n(3)*hh(:,3);

                f = fn .* q / (1j*k0)^2;
            
            % reduced kernels from 1/R     
            case 2
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                
                ff(:,1) = (X-Xp) .* R/8;
                ff(:,2) = (Y-Yp) .* R/8;
                ff(:,3) = (Z-Zp) .* R/8;
                
                q = n(1)*ff(:,1) + n(2)*ff(:,2) + n(3)*ff(:,3);
                
                f = fn .* q / 4/pi;
            
            % reduced kernels from R 
            case 3
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                
                ff(:,1) = (X-Xp) .* R.^3/72;
                ff(:,2) = (Y-Yp) .* R.^3/72;
                ff(:,3) = (Z-Zp) .* R.^3/72;
                
                q = n(1)*ff(:,1) + n(2)*ff(:,2) + n(3)*ff(:,3);
                
                f = fn .* q *(-k0^2/8/pi);
        end
        
    case 3
        %% 3rd surface-surface kernel
        % fm
        switch l
            case 1
                fm = 1;
            case 2
                fm = 1;
            case 3
                fm = X  - r_m(1); 
            case 4
                fm = 1;
            case 5
                fm = Y  - r_m(2);
            case 6
                fm = 1;
            case 7 
                fm = Z  - r_m(3);
            case 8
                fm = X  - r_m(1); 
            case 9
                fm = Y  - r_m(2);
            case 10
                fm = Z  - r_m(3);
        end
        
        switch volume_ker
            % reduced kernels from G
            case 1
                % Green
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                G = 1/4/pi * exp(-1i*k0*R) ./ R;
                G0 = 1/4/pi ./ R;

                % F
                ff(:,1) = (X-Xp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
                ff(:,2) = (Y-Yp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);
                ff(:,3) = (Z-Zp)/(1i*k0)^2 .* ( -1i*k0*G./R - G./R.^2 + G0./R.^2);

                % F - RG0/2
                gg(:,1) = (X-Xp)/2.*G0;
                gg(:,2) = (Y-Yp)/2.*G0;
                gg(:,3) = (Z-Zp)/2.*G0;
                hh = ff - gg;


                q = n(1)*hh(:,1) + n(2)*hh(:,2) + n(3)*hh(:,3);

                f = fm .* q / (1j*k0)^2;
                
            % reduced kernels from 1/R 
            case 2
                 R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                
                ff(:,1) = (X-Xp) .* R/8;
                ff(:,2) = (Y-Yp) .* R/8;
                ff(:,3) = (Z-Zp) .* R/8;
                
                q = n(1)*ff(:,1) + n(2)*ff(:,2) + n(3)*ff(:,3);
                
                f = fm .* q / 4/pi;
            
            % reduced kernels from R 
            case 3
                
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                
                ff(:,1) = (X-Xp) .* R.^3/72;
                ff(:,2) = (Y-Yp) .* R.^3/72;
                ff(:,3) = (Z-Zp) .* R.^3/72;
                
                q = n(1)*ff(:,1) + n(2)*ff(:,2) + n(3)*ff(:,3);
                
                f = fm .* q *(-k0^2/8/pi);
        end
                
        
        
    case 4
        %% 4th surface-surface kernel
        switch volume_ker
            % reduced kernels from G
            case 1
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);

                G = 1/4/pi * exp(-1i*k0*R) ./ R;
                G0 = 1/4/pi./ R;

                f = ( (G-G0)/(1j*k0)^2 - R/8/pi ) / (1j*k0)^2 ;
                
            % reduced kernels from 1/R 
            case 2
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                f = R.^3/24 / 4/pi;
            
            % reduced kernels from R 
            case 3
                R = sqrt((X-Xp).^2+(Y-Yp).^2+(Z-Zp).^2);
                f = R.^5/360 *(-k0^2/8/pi);
        end
        
end
        
        
        
        
end