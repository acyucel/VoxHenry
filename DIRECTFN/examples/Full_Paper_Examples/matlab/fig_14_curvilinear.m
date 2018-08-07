addpath('../../../mex')

ko = 2*pi;


d = 0.1;
dd = d / sqrt(2.0);

Nmax = 25;
Npoints = 1:Nmax;
Error_ST = zeros(Nmax,1);
Error_EA = zeros(Nmax,1);
Error_VA = zeros(Nmax,1);
Nref = 20;

%% Self-Term Case
r1 = [ d, 0.0, 0.0 ];
r2 = [ d, 0.0, 0.5*d ];
r3 = [ d, 0.0, d ];
r4 = [ dd, dd, 0.0 ];
r5 = [ dd, dd, 0.5*d ];
r6 = [ dd, dd, d ];
r7 = [ 0.0, d, 0.0 ];
r8 = [ 0.0, d, 0.5*d ];
r9 = [ 0.0, d, d ];

fprintf('Computing reference value for ST case ...\n');
I_ST_ref = directfn_quad_st_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, Nref, Nref, Nref, Nref, ko, 'Vector_WS');
fprintf('Reference value computed\n');

fprintf('Convergence test starting for ST case...\n');
for N = Npoints
    fprintf('... N = %d\n', N);
    I_ST = directfn_quad_st_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, N, N, N, N, ko, 'Vector_WS');
    Error = abs(abs(I_ST - I_ST_ref)./ abs(I_ST_ref) + eps);
    for i = 1:16
        if abs(I_ST_ref(i)) < eps
            Error(i) = eps;
        end
    end
    Error_ST(N) = max(Error);
end
fprintf('Convergence test completed.\n');
%%
%% Edje Adjacent Case

r1 = [ d, 0.0, 0.0 ];
r2 = [ d, 0.0, 0.5*d ];
r3 = [ d, 0.0, d ];
r4 = [ dd, dd, 0.0 ];
r5 = [ dd, dd, 0.5*d ];
r6 = [ dd, dd, d ];
r7 = [ 0.0, d, 0.0 ];
r8 = [ 0.0, d, 0.5*d ];
r9 = [ 0.0, d, d ];
r10 = [ -dd, dd, 0.0 ];
r11 = [ -dd, dd, 0.5*d ];
r12 = [ -dd, dd, d ];
r13 = [ -d, 0.0, 0.0 ];
r14 = [ -d, 0.0, 0.5*d ];
r15 = [ -d, 0.0, d ];

fprintf('Computing reference value for EA case ...\n');
I_EA_ref = directfn_quad_ea_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, ...
    Nref, Nref, Nref, Nref, ko, 'Vector_SS');
fprintf('Reference value computed\n');

fprintf('Convergence test starting for EA case...\n');
for N = Npoints
    fprintf('... N = %d\n', N);
    I_EA = directfn_quad_ea_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, ...
        N, N, N, N, ko, 'Vector_SS');
    Error = abs(abs(I_EA - I_EA_ref)./ abs(I_EA_ref) + eps);
    for i = 1:16
        if abs(I_EA_ref(i)) < eps
            Error(i) = eps;
        end
    end
    Error_EA(N) = max(Error);
end
fprintf('Convergence test completed.\n');
%%
%% Vertex Adjacent Case
r1 = [d, 0.0, 0.0 ];
r2 = [ d, 0.0, 0.5*d ];
r3 = [ d, 0.0, d ];
r4 = [dd, dd, 0.0 ];
r5 = [ dd, dd, 0.5*d ];
r6 = [ dd, dd, d ];
r7 = [ 0.0, d, 0.0 ];
r8 = [ 0.0, d, 0.5*d ];
r9 = [ 0.0, d, d ];
r10 = [ 0.0, d, 1.5*d ];
r11 = [ 0.0, d, 2.0*d ];
r12 = [ -dd, dd, d ];
r13 = [ -dd, dd, 1.5*d ];
r14 = [ -dd, dd, 2.0*d ];
r15 = [ -d, 0.0, d ];
r16 = [ -d, 0.0, 1.5*d ];
r17 = [ -d, 0.0, 2.0*d ];

fprintf('Computing reference value for VA case ...\n');
I_VA_ref = directfn_quad_va_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, ...
    r16, r17, Nref, Nref, Nref, Nref, ko, 'Vector_SS');
fprintf('Reference value computed\n');

fprintf('Convergence test starting for VA case...\n');
for N = Npoints
    fprintf('... N = %d\n', N);
    I_VA = directfn_quad_va_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, ...
        r16, r17, N, N, N, N, ko, 'Vector_SS');
    Error = abs(abs(I_VA - I_VA_ref)./ abs(I_VA_ref) + eps);
    for i = 1:16
        if abs(I_VA_ref(i)) < eps
            Error(i) = eps;
        end
    end
    Error_VA(N) = max(Error);
end
fprintf('Convergence test completed.\n');
%%
%%
% Plotting Relative Error
figure
semilogy(Npoints, Error_ST, Npoints, Error_EA, ':o', Npoints, Error_VA, '-.s');
grid on;

xlabel('N');
ylabel('relative error');
legend('ST case, weakly singular',' EA case, strongly singular','VA case, strongly singular');