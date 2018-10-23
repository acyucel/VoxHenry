d = 0.1;

ko = 2*pi;

btt = 'Constant';

r1  = [0.0 , 0.0 , 0.0];
r2 = [d , 0.0 , 0.0];




Nmax = 30;
Npoints = 1:Nmax;
M = 8;

% Arrays for errors and timings
Error_ST_quad = zeros(Nmax,M);
Error_ST_tri = zeros(Nmax,M);
%time_ST_quad = zeros(Nmax,1);
%time_ST_tri = zeros(Nmax,1);

%Error_EA_quad = zeros(Nmax,1);
%Error_EA_tri = zeros(Nmax,1);
%time_EA_quad = zeros(Nmax,1);
%time_EA_tri = zeros(Nmax,1);

%Error_VA_quad = zeros(Nmax,1);
%Error_VA_tri = zeros(Nmax,1);
%time_VA_quad = zeros(Nmax,1);
%time_VA_tri = zeros(Nmax,1);
step = 0.2290;
q = zeros(M,1);

for m = 1:M
    
    fprintf('m = %d\n',m);
    
    r3 = [d + step*d*(m-1),  d , 0.0];
    r4 = [0.0+ step*d*(m-1) , d , 0.0];
    q(m,1) = atan(step*(m-1))*2/pi;
    fprintf('Computing reference values ...\n');
    I_ST_ref_quad = directfn_quad_st_plan(r1, r2, r3, r4, Nmax,Nmax,Nmax,Nmax, ko, btt);
    
    %I1 = directfn_tri_st_plan(r1, r2, r4, 25, 25, 25, 25, ko, btt);
    %I2 = directfn_tri_ea_plan(r4, r2, r3, r1, 25, 25, 25, 25, ko, btt);
    %I3 = directfn_tri_ea_plan(r2, r4, r1, r3, 25, 25, 25, 25, ko, btt);
    %I4 = directfn_tri_st_plan(r4, r2, r3, 25, 25, 25, 25, ko, btt);
    %I_ST_ref_tri = I1(1) + I2(1) + I3(1) + I4(1);
    fprintf('Reference values computed.\n');
    for N = Npoints
    
        I_ST_quad = directfn_quad_st_plan(r1, r2, r3, r4, N, N, N, N, ko, btt);
    
        Error_ST_quad(N,m) = abs(abs((I_ST_quad(1) - I_ST_ref_quad(1))) / abs(I_ST_ref_quad(1)) + eps);
    
    
        % directfn_tri
    
        %I1 = directfn_tri_st_plan(r1, r2, r4, N, N, N, N, ko, btt);
        %I2 = directfn_tri_ea_plan(r4, r2, r3, r1, N, N, N, N, ko, btt);
        %I3 = directfn_tri_ea_plan(r2, r4, r1, r3, N, N, N, N, ko, btt);
        %I4 = directfn_tri_st_plan(r4, r2, r3, N, N, N, N, ko, btt);
        %I_ST_tri = I1(1) + I2(1) + I3(1) + I4(1);
        %Error_ST_tri(N,m) = abs(abs((I_ST_tri - I_ST_ref_tri(1))) / abs(I_ST_ref_tri(1)) + eps);
    
    end
end

% Plotting Relative Error


Q = strcat('s = ',num2str(q));
figure(1);
for m = 1:M
    semilogy(Npoints, Error_ST_quad(:,m)); hold on;
    legend(Q);
end 

%figure(2);
%for m = 1:M
    %semilogy(Npoints, Error_ST_tri(:,m), styles(m)); hold on;
%end

ylim([10^-16 10^0]);
grid on;
xlabel('N');
ylabel('relative error');







