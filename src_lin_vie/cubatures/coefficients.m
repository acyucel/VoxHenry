function f = coefficients(n,np,l,k0,dx,ker_type)

%% calculates the constnat coefficients part for the 4 surface-surface kernels
switch ker_type
    
    case 1
        %% 1st kernel's coefficients
        f = -dot(n,np);

    case 2
        %% 2nd kernel's coefficients
        % gradinet of fm
        switch l
            case 1
                Nfm = [0 0 0];
            case 2
                Nfm = [0 0 0];
            case 3
                Nfm = [1 0 0];
            case 4
                Nfm = [0 0 0];
            case 5
                Nfm = [0 1 0];
            case 6
                Nfm = [0 0 0];
            case 7 
                Nfm = [0 0 1];
            case 8
                Nfm = [1 0 0];
            case 9
                Nfm = [0 1 0];
            case 10
                Nfm = [0 0 1];
        end

        f = dot(np,Nfm);
        
    case 3
        %% 3rd kernel's coefficients
        % gradient of fn
        switch l
            case 1
                Nfn = [0 0 0];
            case 2
                Nfn = [1 0 0];
            case 3
                Nfn = [0 0 0];
            case 4
                Nfn = [0 1 0];
            case 5
                Nfn = [0 0 0];
            case 6
                Nfn = [0 0 1];
            case 7 
                Nfn = [0 0 0];
            case 8
                Nfn = [1 0 0];
            case 9
                Nfn = [0 1 0];
            case 10
                Nfn = [0 0 1];
        end


        f = -dot(np,Nfn);

        
    case 4
        %% 4th kernel's coefficients
        switch l
    
            case 1
                Nfm = [0 0 0];
                Nfn = [0 0 0];
            case 2
                Nfm = [0 0 0];
                Nfn = [1 0 0];
            case 3
                Nfm = [1 0 0];
                Nfn = [0 0 0];
            case 4
                Nfm = [0 0 0];
                Nfn = [0 1 0];
            case 5
                Nfm = [0 1 0];
                Nfn = [0 0 0];
            case 6
                Nfm = [0 0 0];
                Nfn = [0 0 1];
            case 7 
                Nfm = [0 0 1];
                Nfn = [0 0 0];
            case 8
                Nfm = [1 0 0];
                Nfn = [1 0 0];
            case 9
                Nfm = [0 1 0];
                Nfn = [0 1 0];
            case 10
                Nfm = [0 0 1];
                Nfn = [0 0 1];

        end

        f = dot(np,Nfn) * dot(n,Nfm);
        
end


end
