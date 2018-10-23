addpath('../../../mex')

d = 0.1;

ko = 2*pi;

btt = 'Constant';

r1 = [ 0.0 , 0.4*d , 0.0 ];
r2 = [ d , 0.0 , 0.0 ];
r3 = [ 0.8*d,  d , 0.0 ];
r4 = [ 0.2*d , 0.8*d , 0.0 ];

fprintf('Computing reference value ...\n');
I_ST_ref = directfn_quad_st_plan(r1, r2, r3, r4, 20, 20, 20, 20, ko, btt);
fprintf('Reference value computed\n');

Nmax = 16;
Npoints = 1:Nmax;

Error = zeros(Nmax, 1);
time = zeros(Nmax, 1);

Error_1 = zeros(Nmax, 1);
time_1 = zeros(Nmax, 1);

Error_2 = zeros(Nmax, 1);
time_2 = zeros(Nmax, 1);

Error_3 = zeros(Nmax, 1);
time_3 = zeros(Nmax, 1);

Error_4 = zeros(Nmax, 1);
time_4 = zeros(Nmax, 1);


for N = Npoints
    fprintf('... N = %d\n', N);
    % Same N
    tic;
    I_ST_quad = directfn_quad_st_plan(r1, r2, r3, r4, N, N, N, N, ko, btt);
    Error(N) = abs(abs((I_ST_quad(1) - I_ST_ref(1))) / abs(I_ST_ref(1)) + eps);
    time(N) = toc;
    
    
    % Varying N1
    tic;
    I_ST_quad = directfn_quad_st_plan(r1, r2, r3, r4, N, 20, 20, 20, ko, btt);
    Error_1(N) = abs(abs((I_ST_quad(1) - I_ST_ref(1))) / abs(I_ST_ref(1)) + eps);
    time_1(N) = toc;
    
    % Varying N2
    tic;
    I_ST_quad = directfn_quad_st_plan(r1, r2, r3, r4, 20, N, 20, 20, ko, btt);
    Error_2(N) = abs(abs((I_ST_quad(1) - I_ST_ref(1))) / abs(I_ST_ref(1)) + eps);
    time_2(N) = toc;
    
    % Varying N3
    tic;
    I_ST_quad = directfn_quad_st_plan(r1, r2, r3, r4, 20, 20, N, 20, ko, btt);
    Error_3(N) = abs(abs((I_ST_quad(1) - I_ST_ref(1))) / abs(I_ST_ref(1)) + eps);
    time_3(N) = toc;
    
    % Varying N4
    tic;
    I_ST_quad = directfn_quad_st_plan(r1, r2, r3, r4, 20, 20, 20, N, ko, btt);
    Error_4(N) = abs(abs((I_ST_quad(1) - I_ST_ref(1))) / abs(I_ST_ref(1)) + eps);
    time_4(N) = toc;
end



% Plotting Relstive Error
figure
semilogy(Npoints, Error, Npoints, Error_1, ':o', Npoints, Error_2, '-.s', Npoints, Error_3, ':v', Npoints, Error_4, '-.s');
grid on;

xlabel('N');
ylabel('relative error');
legend('Same N','Varying N_1', 'Varying N_2', 'Varying N_3', 'Varying N_4');





