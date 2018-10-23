d = 0.1;

ko = 2*pi;

btt = 'Constant';

r1  = [0.0 , 0.0 , 0.0];
r2 = [d , 0.0 , 0.0];
r3 = [d,  d , 0.0];
r4 = [0.0 , d , 0.0];
r5 = [0.0 , 2.0*d , 0.0];
r6 = [d , 2.0*d , 0.0];
r7 = [2.0*d , 2.0*d , 0.0];
r8 = [2.0*d , d, 0.0];
fprintf('Computing reference values ...\n');
I_ST_ref = directfn_quad_st_plan(r1, r2, r3, r4, 25, 25, 25, 25, ko, btt);
I_EA_ref = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, 25, 25, 25, 25, ko, btt);
I_VA_ref = directfn_quad_va_plan(r1, r2, r3, r4, r8, r7, r6, 25, 25, 25, 25, ko, btt);
fprintf('Reference values computed.\n');

Nmax = 25;
Npoints = 1:Nmax;

% Arrays for errors and timings
Error_ST_quad = zeros(Nmax,1);
Error_ST_tri = zeros(Nmax,1);
time_ST_quad = zeros(Nmax,1);
time_ST_tri = zeros(Nmax,1);

Error_EA_quad = zeros(Nmax,1);
Error_EA_tri = zeros(Nmax,1);
time_EA_quad = zeros(Nmax,1);
time_EA_tri = zeros(Nmax,1);

Error_VA_quad = zeros(Nmax,1);
Error_VA_tri = zeros(Nmax,1);
time_VA_quad = zeros(Nmax,1);
time_VA_tri = zeros(Nmax,1);

for N = Npoints
    if N == 1
        Ntimes = 100;
    else
        if  N < 5   
            Ntimes = 10;
        else
            if N >= 5 
                Ntimes = 3;
            end
        end
    end
    fprintf('... N = %d\n', N);
    %------ Self-Term Case--------------%
    % directfn_quad
    fprintf('test number = ');
    tic;
    for i = 1:Ntimes
         fprintf('%d ', i);
    I_ST_quad = directfn_quad_st_plan(r1, r2, r3, r4, N, N, N, N, ko, btt);
    end
    time_ST_quad(N) = toc/Ntimes;
    fprintf('\n');
    Error_ST_quad(N) = abs(abs((I_ST_quad(1) - I_ST_ref(1))) / abs(I_ST_ref(1)) + eps);
    
    
    % directfn_tri
    tic;
    I1 = directfn_tri_st_plan(r1, r2, r4, N, N, N, N, ko, btt);
	I2 = directfn_tri_ea_plan(r4, r2, r3, r1, N, N, N, N, ko, btt);
	I3 = directfn_tri_ea_plan(r2, r4, r1, r3, N, N, N, N, ko, btt);
	I4 = directfn_tri_st_plan(r4, r2, r3, N, N, N, N, ko, btt);
    I_ST_tri = I1(1) + I2(1) + I3(1) + I4(1);
    Error_ST_tri(N) = abs(abs((I_ST_tri - I_ST_ref(1))) / abs(I_ST_ref(1)) + eps);
    time_ST_tri(N) = toc;
    
    %------ Edge Adjacent Case--------------%
    
    % directfn_quad
    tic;
    I_EA_quad = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, N, N, N, N, ko, btt);
    Error_EA_quad(N) = abs(abs((I_EA_quad(1) - I_EA_ref(1))) / abs(I_EA_ref(1)) + eps);
    time_EA_quad(N) = toc;
    
    % directfn_tri
    tic;
    I1 = directfn_tri_ea_plan(r4, r3, r6, r2, N, N, N, N, ko, btt);
	I2 = directfn_tri_va_plan(r4, r2, r3, r5, r6, N, N, N, N, ko, btt);
	I3 = directfn_tri_va_plan(r4, r1, r2, r5, r6, N, N, N, N, ko, btt);
	I4 = directfn_tri_va_plan(r4, r1, r2, r6, r3, N, N, N, N, ko, btt);
    
    I_EA_tri = I1(1) + I2(1) + I3(1) + I4(1);
    Error_EA_tri(N) = abs(abs((I_EA_tri - I_EA_ref(1))) / abs(I_EA_ref(1)) + eps);
    time_EA_tri(N) = toc;
    
    
    %------ Vertex Adjacent Case--------------%
    
    % directfn_quad
    tic;
    I_VA_quad = directfn_quad_va_plan(r1, r2, r3, r4, r8, r7, r6, N, N, N, N, ko, btt);
    Error_VA_quad(N) = abs(abs((I_VA_quad(1) - I_VA_ref(1))) / abs(I_VA_ref(1)) + eps);
    time_VA_quad(N) = toc;
    
    % directfn_tri
    tic;
    I1 = directfn_tri_va_plan(r3, r7, r6, r2, r1, N, N, N, N, ko, btt);
	I2 = directfn_tri_va_plan(r3, r7, r6, r1, r4, N, N, N, N, ko, btt);
	I3 = directfn_tri_va_plan(r3, r8, r7, r1, r4, N, N, N, N, ko, btt);
	I4 = directfn_tri_va_plan(r3, r8, r7, r2, r1, N, N, N, N, ko, btt);
    
    I_VA_tri = I1(1) + I2(1) + I3(1) + I4(1);
    Error_VA_tri(N) = abs(abs((I_VA_tri - I_VA_ref(1))) / abs(I_VA_ref(1)) + eps);
    time_VA_tri(N) = toc;
end

% Plotting Relative Error
figure;
semilogy(Npoints, Error_ST_quad, '--o', Npoints, Error_ST_tri, ':o'); hold on;
semilogy(Npoints, Error_EA_quad, '--v', Npoints, Error_EA_tri, ':v'); hold on;
semilogy(Npoints, Error_VA_quad, '--*', Npoints, Error_VA_tri, ':*'); hold on;
grid on;
title('Quadrilateral VS Triangular');
xlabel('N');
ylabel('relative error');
lbl1 = 'DIRECTFN-quad, ST case';
lbl2 = 'DIRECTFN-tri, ST case';
lbl3 = 'DIRECTFN-quad, EA case';
lbl4 = 'DIRECTFN-tri, EA case';
lbl5 = 'DIRECTFN-quad, VA case';
lbl6 = 'DIRECTFN-tri, VA case';
legend(lbl1, lbl2, lbl3, lbl4, lbl5, lbl6);

% Plotting CPU time
figure
plot(Npoints, time_ST_tri./time_ST_quad, '-o', Npoints, time_EA_tri./time_EA_quad, '-s', Npoints, time_VA_tri./time_VA_quad, '-^');
grid on;
title('CPU time');
xlabel('N');
ylabel('$$\frac{T_{tri}}{T_{quad}}$$','Interpreter','Latex');
legend('ST case','EA case','VA case');





