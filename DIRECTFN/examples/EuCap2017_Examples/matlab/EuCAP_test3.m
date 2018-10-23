addpath('../../../mex')

ko = 2*pi;


d = 0.1;
dd = d / sqrt(2.0);

Nmax = 25;
Npoints = 1:Nmax;
Nref = 20;
Error = zeros(Nmax,1);

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
I_ST_ref = directfn_quad_st_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, Nref, Nref, Nref, Nref, ko, 'Constant');
fprintf('Reference value computed\n');

fprintf('Convergence test starting for ST case...\n');
for N = Npoints
    fprintf('... N = %d\n', N);
    I_ST = directfn_quad_st_curv(r1, r2, r3, r4, r5, r6, r7, r8, r9, N, N, N, N, ko, 'Constant');
    Error(N) = abs(abs(I_ST - I_ST_ref)/ abs(I_ST_ref) + eps);
end
fprintf('Convergence test completed.\n');

% Plotting Relative Error
figure
semilogy(Npoints, Error);
grid on;

xlabel('N');
ylabel('relative error');
