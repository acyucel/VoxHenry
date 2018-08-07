d = 0.1;

ko = 2*pi;

btt = 'Vector_WS';

r1 = [ 0.0 , 0.0 , 0.0];
r2 = [ d , 0.0 , 0.0];
r3 = [ d,  d , 0.0];
r4 = [0.0 , d , 0.0];
r5 = [ 2.0*d, 0.0 , 0.0];
r6 = [ 2.0*d ,d , 0.0];
r7 = [ d , 0.0 , d];
r8 = [ d , d , d];
r9 = [d , 2.0*d , d];
r10 = [ d , 2.0*d , 0.0];
r11 = [ 2.0* d , 2.0*d , 0.0];

Nref = 20;
fprintf('Computing reference value ...\n');
I_ST_ref = directfn_quad_st_plan (r1, r2, r3, r4, Nref, Nref, Nref, Nref, ko, btt);

fprintf('Reference value computed\n');

Nmax = 20;
Npoints = 1:Nmax;

Error_0 = zeros(Nmax,1);

Error_1 = zeros(Nmax,1);

Error_2 = zeros(Nmax,1);

Error_3 = zeros(Nmax,1);

Error_4 = zeros(Nmax,1);


for N = Npoints
    fprintf('... N = %d\n', N);
    
    % Same N
    fprintf('...N N N N\n');
    I_ST = directfn_quad_st_plan(r1, r2, r3, r4, N, N, N, N, ko, btt);
    Error_0(N) = max(abs(abs(I_ST - I_ST_ref)./ abs(I_ST_ref) + eps));
    
    
    % N1
    fprintf('...N %d %d %d\n', Nref, Nref, Nref);
    I_ST = directfn_quad_st_plan(r1, r2, r3, r4, N, 20, 20, 20, ko, btt);
    Error_1(N) = max(abs(abs(I_ST - I_ST_ref)./ abs(I_ST_ref) + eps));
    
    % N2
    fprintf('.. %d N %d %d\n', Nref, Nref, Nref);
    I_ST = directfn_quad_st_plan(r1, r2, r3, r4, 20, N, 20, 20, ko, btt);
    Error_2(N) = max(abs(abs(I_ST - I_ST_ref)./ abs(I_ST_ref) + eps));
    
    
    % N3
    fprintf('...%d %d N %d\n', Nref, Nref, Nref);
    I_ST = directfn_quad_st_plan(r1, r2, r3, r4, 20, 20, N, 20, ko, btt);
    Error_3(N) = max(abs(abs(I_ST - I_ST_ref)./ abs(I_ST_ref) + eps));
    
    
    % N4
    fprintf('...%d %d %d N\n', Nref, Nref, Nref);
    I_ST = directfn_quad_st_plan(r1, r2, r3, r4, 20, 20, 20, N, ko, btt);
    Error_4(N) = max(abs(abs(I_ST - I_ST_ref)./ abs(I_ST_ref) + eps));
    
end



% Plotting Relstive Error
figure
semilogy(Npoints, Error_0, Npoints, Error_1, ':o', Npoints, Error_2, '-.s', Npoints, Error_3, ':v', Npoints, Error_4, '-.s');
grid on;

xlabel('N');
ylabel('relative error');
legend('Same N','Varying N_1', 'Varying N_2', 'Varying N_3', 'Varying N_4');




