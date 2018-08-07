addpath('../../../mex')
% add path or so (with mexed binary located at ROOT_DIR/mex/)

Ntimes = 10;

d = 0.1;

ko = 2*pi;

btt = 'Constant';

r1 = [ 0.0 , 0.4*d , 0.0 ];
r2 = [ d , 0.0 , 0.0 ];
r3 = [ 0.8*d,  d , 0.0 ];
r4 = [ 0.2*d , 0.8*d , 0.0 ];

fprintf('Computing reference value ...\n');
I_ST_ref = directfn_quad_st_plan(r1, r2, r3, r4, 25, 25, 25, 25, ko, btt);
fprintf('Reference value computed\n');

Nmax = 30;
Npoints = 1:Nmax;

Error_quad = zeros(Nmax,1);
Error_tri = zeros(Nmax,1);
time_quad = zeros(Nmax,1);
time_tri = zeros(Nmax,1);
for N = Npoints
    fprintf('... N = %d\n', N);
    % directfn_quad
    tic;
    for i = 1:Ntimes
        I_ST_quad = directfn_quad_st_plan(r1, r2, r3, r4, N, N, N, N, ko, btt);
    end
    time_quad(N) = toc/Ntimes;
    Error_quad(N) = abs(abs((I_ST_quad(1) - I_ST_ref(1))) / abs(I_ST_ref(1)) + eps);
    
    
    % directfn_tri
    tic;
    for i = 1:Ntimes
        I1 = directfn_tri_st_plan(r1, r2, r4, N, N, N, N, ko, btt);
        I2 = directfn_tri_ea_plan(r4, r2, r3, r1, N, N, N, N, ko, btt);
        I3 = directfn_tri_ea_plan(r2, r4, r1, r3, N, N, N, N, ko, btt);
        I4 = directfn_tri_st_plan(r4, r2, r3, N, N, N, N, ko, btt);
        I_ST_tri = I1(1) + I2(1) + I3(1) + I4(1);
    end
    time_tri(N) = toc/Ntimes;
    Error_tri(N) = abs(abs((I_ST_tri - I_ST_ref(1))) / abs(I_ST_ref(1)) + eps);
    
    
end

% Plotting Relstive Error
figure;
semilogy(Npoints, Error_quad, ':o', Npoints, Error_tri, '-.s');
grid on;
title('Quadrilateral VS Triangular');
xlabel('N');
ylabel('relative error');
legend('DIRECTFN-quad','DIRECTFN-tri');

% Plotting CPU time
figure('position',[200, 100, 700, 350]);
subplot(121);
plot(Npoints, time_quad, ':o', Npoints, time_tri, '-.s');
grid on;
title('CPU time');
xlabel('N');
ylabel('time');
legend('DIRECTFN-quad','DIRECTFN-tri');

subplot(122);
reltime = time_tri./time_quad;
plot(Npoints, reltime, '-o');
grid on;
title('CPU time');
xlabel('N');
str = '$$\frac{t_{tri}}{t_{quad}}$$';
ylabel(str,'Interpreter','Latex');



