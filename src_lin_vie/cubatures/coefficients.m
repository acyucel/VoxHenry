function f = coefficients(n,np,l,k0,dx,ker_type)

%% calculates the constant coefficients part for the 4 surface-surface kernels
% for all the cases of Table 1,2,3 the coefficients are the same.
% they differ along the rows all tables

% n: normal vector to each face of the observation voxel 
% np: normal vector to each face of the source voxel
% l: indexing of coefficients according to equations 1a-1j of the report
% ker_type: indexing of the surface-surface kernels according to the rows of tables
% 1,2,3 of the report


switch ker_type
    
    case 1
        %% 1st kernel's coefficients (first row of the tables)
        f = -dot(n,np);

    case 2
        %% 2nd kernel's coefficients (second row of the tables)
        % gradient of fm
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
        %% 3rd kernel's coefficients (third row of the tables)
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
        %% 4th kernel's coefficients (fourth row of the tables)
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
