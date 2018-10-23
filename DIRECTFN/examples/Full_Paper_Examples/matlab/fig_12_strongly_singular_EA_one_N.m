d = 0.1;

ko = 2*pi;

btt = 'Vector_SS';

r1 = [ 0.0 , 0.0 , 0.0];
r2 = [ d , 0.0 , 0.0 ];
r3 = [ d,  d , 0.0 ];
r4 = [ 0.0 , d , 0.0 ];
r5 = [ d , 0.0 , d ];
r6 = [ d , d , d ];

Nref = 20;
fprintf('Computing reference value ...\n');
I_EAo_ref = directfn_quad_ea_plan ( r1, r2, r3, r4, r5, r6, Nref, Nref, Nref, Nref, ko, btt);


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
    I_EAo = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, N, N, N, N, ko, btt);
    Error = abs(abs(I_EAo - I_EAo_ref)./ abs(I_EAo_ref) + eps);
    for i = 1:16
        if abs(I_EAo_ref(i)) < eps
            Error(i) = eps;
        end
    end
    Error_0(N) = max(Error);
    
    
    % N1
    fprintf('...N %d %d %d\n', Nref, Nref, Nref);
    I_EAo = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, N, 20, 20, 20, ko, btt);
    Error = abs(abs(I_EAo - I_EAo_ref)./ abs(I_EAo_ref) + eps);
    for i = 1:16
        if abs(I_EAo_ref(i)) < eps
            Error(i) = eps;
        end
    end
    Error_1(N) = max(Error);
    
    % N2
    fprintf('.. %d N %d %d\n', Nref, Nref, Nref);
    I_EAo = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, 20, N, 20, 20, ko, btt);
    Error = abs(abs(I_EAo - I_EAo_ref)./ abs(I_EAo_ref) + eps);
    for i = 1:16
        if abs(I_EAo_ref(i)) < eps
            Error(i) = eps;
        end
    end
    Error_2(N) = max(Error);
    
    
    % N3
    fprintf('...%d %d N %d\n', Nref, Nref, Nref);
    I_EAo = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, 20, 20, N, 20, ko, btt);
    Error = abs(abs(I_EAo - I_EAo_ref)./ abs(I_EAo_ref) + eps);
    for i = 1:16
        if abs(I_EAo_ref(i)) < eps
            Error(i) = eps;
        end
    end
    Error_3(N) = max(Error);
    
    
    % N4
    fprintf('...%d %d %d N\n', Nref, Nref, Nref);
    I_EAo = directfn_quad_ea_plan(r1, r2, r3, r4, r5, r6, 20, 20, 20, N, ko, btt);
    Error = abs(abs(I_EAo - I_EAo_ref)./ abs(I_EAo_ref) + eps);
    for i = 1:16
        if abs(I_EAo_ref(i)) < eps
            Error(i) = eps;
        end
    end
    Error_4(N) = max(Error);
    
end



% Plotting Relative Error
figure
semilogy(Npoints, Error_0, Npoints, Error_1, ':o', Npoints, Error_2, '-.s', Npoints, Error_3, ':v', Npoints, Error_4, '-.s');
grid on;

xlabel('N');
ylabel('relative error');
legend('Same N','Varying N_1', 'Varying N_2', 'Varying N_3', 'Varying N_4');




