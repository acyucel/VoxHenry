function [adjacency_type,ordered_points] = points_mappping(m,dx)

    % adjacency_type: [6x6]     matrix with the singularity type 
    % orderd_points : [3x7x6x6] matrix with the 7 3D points for every m and
    % adjacency type. For ST the last 3 points are zero, for EA the last 2.


    % cases for the singularity type for surface-surface integrals over squares
    % ST  -> self term
    % EAo -> edge adjacent orthogonal
    % EAc -> edge adjacent coplanar
    % VAo -> vertex adjacent orthogonal
    % VAc -> vertex adjacent coplanar
    % 0.0 -> non-singular integral
    ST  = 1;
    EAo = 2;
    EAc = 2;
    VAo = 3;
    VAc = 3;

    a = -dx/2;
    b = +dx/2;
    % primed points - hold stable - basis voxel withe centre (0,0,0)
    p1 = [a a a]';
    p2 = [b a a]';
    p3 = [b b a]';
    p4 = [a b a]';
    p5 = [a a b]';
    p6 = [b a b]';
    p7 = [b b b]';
    p8 = [a b b]';

    % unprimed points - moved - testing voxel with centre dr
    d = [dx dx dx];
    % distance vector
    dr = ((m-1) .* d)';
    r1 = p1 + dr;
    r2 = p2 + dr;
    r3 = p3 + dr;
    r4 = p4 + dr;
    r5 = p5 + dr;
    r6 = p6 + dr;
    r7 = p7 + dr;
    r8 = p8 + dr;

    % zero vector
    ze = zeros(3,1);
    % n  = [-x,x,-y,y,-z,z]
    % n' = [-x,x,-y,y,-z,z]

    % ({x,y,z},r1...r7,n,n')
    ordered_points = NaN(3,7,6,6);
    % (n,n')
    adjacency_type = zeros(6,6);

    if isequal(m,[1,1,1]) % m == [1,1,1]

        ordered_points(:,:,1,1) = [p1 p4 p8 p5 ze ze ze];
        ordered_points(:,:,1,3) = [p2 p6 p5 p1 p4 p8 ze];
        ordered_points(:,:,1,4) = [p3 p7 p8 p4 p1 p5 ze];
        ordered_points(:,:,1,5) = [p2 p3 p4 p1 p5 p8 ze];
        ordered_points(:,:,1,6) = [p6 p7 p8 p5 p1 p4 ze];

        ordered_points(:,:,2,2) = [p2 p3 p7 p6 ze ze ze];
        ordered_points(:,:,2,3) = [p1 p5 p6 p2 p3 p7 ze];
        ordered_points(:,:,2,4) = [p4 p8 p7 p3 p2 p6 ze];
        ordered_points(:,:,2,5) = [p1 p4 p3 p2 p6 p7 ze];
        ordered_points(:,:,2,6) = [p8 p5 p6 p7 p3 p2 ze];

        ordered_points(:,:,3,1) = [p4 p8 p5 p1 p2 p6 ze];
        ordered_points(:,:,3,2) = [p3 p7 p6 p2 p1 p5 ze];
        ordered_points(:,:,3,3) = [p1 p2 p6 p5 ze ze ze];
        ordered_points(:,:,3,5) = [p3 p4 p1 p2 p6 p5 ze];
        ordered_points(:,:,3,6) = [p7 p8 p5 p6 p2 p1 ze];

        ordered_points(:,:,4,1) = [p1 p5 p8 p4 p3 p7 ze];
        ordered_points(:,:,4,2) = [p2 p6 p7 p3 p4 p8 ze];
        ordered_points(:,:,4,4) = [p4 p3 p7 p8 ze ze ze];
        ordered_points(:,:,4,5) = [p1 p2 p3 p4 p8 p7 ze];
        ordered_points(:,:,4,6) = [p5 p6 p7 p8 p4 p3 ze];

        ordered_points(:,:,5,1) = [p8 p5 p1 p4 p3 p2 ze];
        ordered_points(:,:,5,2) = [p6 p7 p3 p2 p1 p4 ze];
        ordered_points(:,:,5,3) = [p5 p6 p2 p1 p4 p3 ze];
        ordered_points(:,:,5,4) = [p7 p8 p4 p3 p2 p1 ze];
        ordered_points(:,:,5,5) = [p1 p2 p3 p4 ze ze ze];

        ordered_points(:,:,6,1) = [p1 p4 p8 p5 p6 p7 ze];
        ordered_points(:,:,6,2) = [p2 p3 p7 p6 p5 p8 ze];
        ordered_points(:,:,6,3) = [p1 p2 p6 p5 p8 p7 ze];
        ordered_points(:,:,6,4) = [p4 p3 p7 p8 p5 p6 ze];
        ordered_points(:,:,6,6) = [p5 p6 p7 p8 ze ze ze];


         adjacency_type = [ ST  0.0 EAo EAo EAo EAo
                            0.0 ST  EAo EAo EAo EAo
                            EAo EAo ST  0.0 EAo EAo
                            EAo EAo 0.0 ST  EAo EAo
                            EAo EAo EAo EAo ST  0.0
                            EAo EAo EAo EAo 0.0 ST ];

    elseif isequal(m,[2,1,1]) % m == [2,1,1]

        ordered_points(:,:,1,2) = [p2 p3 p7 p6 ze ze ze];
        ordered_points(:,:,1,3) = [p1 p5 p6 p2 p3 p7 ze];
        ordered_points(:,:,1,4) = [p4 p8 p7 p3 p2 p6 ze];
        ordered_points(:,:,1,5) = [p1 p4 p3 p2 p6 p7 ze];
        ordered_points(:,:,1,6) = [p5 p8 p7 p6 p2 p3 ze];

        ordered_points(:,:,3,2) = [p3 p7 p6 p2 r2 r6 ze];
        ordered_points(:,:,3,3) = [p1 p5 p6 p2 r2 r6 ze];
        ordered_points(:,:,3,5) = [p4 p1 p2 p3 r5 r6 r2];
        ordered_points(:,:,3,6) = [p8 p5 p6 p7 r1 r2 r6];

        ordered_points(:,:,4,2) = [p2 p6 p7 p3 r3 r7 ze];
        ordered_points(:,:,4,4) = [p4 p8 p7 p3 r3 r7 ze];
        ordered_points(:,:,4,5) = [p1 p4 p3 p2 r8 r7 r3];
        ordered_points(:,:,4,6) = [p5 p8 p7 p6 r4 r3 r7];

        ordered_points(:,:,5,2) = [p6 p7 p3 p2 r2 r3 ze];
        ordered_points(:,:,5,3) = [p5 p6 p2 p1 r2 r3 r4];
        ordered_points(:,:,5,4) = [p8 p7 p3 p4 r3 r2 r1];
        ordered_points(:,:,5,5) = [p1 p4 p3 p2 r2 r3 ze];

        ordered_points(:,:,6,2) = [p2 p3 p7 p6 r6 r7 ze];
        ordered_points(:,:,6,3) = [p1 p5 p6 p2 r8 r7 r6];
        ordered_points(:,:,6,4) = [p4 p8 p7 p3 r5 r6 r7];
        ordered_points(:,:,6,6) = [p5 p8 p7 p6 r6 r7 ze];


        adjacency_type =  [ 0.0 ST  EAo EAo EAo EAo
                            0.0 0.0 0.0 0.0 0.0 0.0
                            0.0 EAo EAc 0.0 VAo VAo
                            0.0 EAo 0.0 EAc VAo VAo
                            0.0 EAo VAo VAo EAc 0.0
                            0.0 EAo VAo VAo 0.0 EAc ];

     elseif isequal(m,[2,1,2]) % m == [2,1,2] 

        ordered_points(:,:,1,2) = [p2 p3 p7 p6 r5 r8 ze];
        ordered_points(:,:,1,3) = [p1 p5 p6 p2 r5 r8 r4];
        ordered_points(:,:,1,4) = [p4 p8 p7 p3 r8 r5 r1];
        ordered_points(:,:,1,6) = [p5 p8 p7 p6 r5 r8 ze];

        ordered_points(:,:,3,2) = [p3 p7 p6 p2 r5 r6 r2];
        ordered_points(:,:,3,3) = [p1 p2 p6 p5 r2 r6 r5];
        ordered_points(:,:,3,6) = [p8 p7 p6 p5 r2 r6 r5];

        ordered_points(:,:,4,2) = [p2 p6 p7 p3 r8 r7 r3];
        ordered_points(:,:,4,4) = [p4 p3 p7 p8 r3 r7 r8];
        ordered_points(:,:,4,6) = [p5 p6 p7 p8 r3 r7 r8];

        ordered_points(:,:,5,2) = [p2 p3 p7 p6 r2 r3 ze];
        ordered_points(:,:,5,3) = [p1 p2 p6 p5 r2 r3 r4];
        ordered_points(:,:,5,4) = [p4 p3 p7 p8 r3 r2 r1];
        ordered_points(:,:,5,6) = [p5 p8 p7 p6 r2 r3 ze];


        adjacency_type =  [ 0.0 EAc VAo VAo 0.0 EAo
                            0.0 0.0 0.0 0.0 0.0 0.0
                            0.0 VAo VAc 0.0 0.0 VAo
                            0.0 VAo 0.0 VAc 0.0 VAo
                            0.0 EAo VAo VAo 0.0 EAc
                            0.0 0.0 0.0 0.0 0.0 0.0 ];

     elseif isequal(m,[1,1,2]) % m == [1,1,2] 

         ordered_points(:,:,1,1) = [p1 p4 p8 p5 r5 r8 ze];
         ordered_points(:,:,1,3) = [p2 p6 p5 p1 r5 r8 r4];
         ordered_points(:,:,1,4) = [p3 p7 p8 p4 r8 r5 r1];
         ordered_points(:,:,1,6) = [p6 p7 p8 p5 r5 r8 ze];

         ordered_points(:,:,2,2) = [p2 p3 p7 p6 r6 r7 ze];
         ordered_points(:,:,2,3) = [p1 p5 p6 p2 r6 r7 r3];
         ordered_points(:,:,2,4) = [p4 p8 p7 p3 r7 r6 r2];
         ordered_points(:,:,2,6) = [p5 p8 p7 p6 r6 r7 ze];

         ordered_points(:,:,3,1) = [p4 p8 p5 p1 r5 r6 r2];
         ordered_points(:,:,3,2) = [p3 p7 p6 p2 r6 r5 r1];
         ordered_points(:,:,3,3) = [p1 p2 p6 p5 r5 r6 ze];
         ordered_points(:,:,3,6) = [p8 p7 p6 p5 r5 r6 ze];

         ordered_points(:,:,4,1) = [p1 p5 p8 p4 r8 r7 r3];
         ordered_points(:,:,4,2) = [p2 p6 p7 p3 r7 r8 r4];
         ordered_points(:,:,4,4) = [p4 p3 p7 p8 r8 r7 ze];
         ordered_points(:,:,4,6) = [p5 p6 p7 p8 r8 r7 ze];

         ordered_points(:,:,5,1) = [p1 p4 p8 p5 p6 p7 ze];
         ordered_points(:,:,5,2) = [p2 p3 p7 p6 p5 p8 ze];
         ordered_points(:,:,5,3) = [p1 p2 p6 p5 p8 p7 ze];
         ordered_points(:,:,5,4) = [p4 p3 p7 p8 p5 p6 ze];
         ordered_points(:,:,5,6) = [p5 p6 p7 p8 ze ze ze];


         adjacency_type = [ EAc 0.0 VAo VAo 0.0 EAo
                            0.0 EAc VAo VAo 0.0 EAo
                            VAo VAo EAc 0.0 0.0 EAo
                            VAo VAo 0.0 EAc 0.0 EAo
                            EAo EAo EAo EAo 0.0 ST
                            0.0 0.0 0.0 0.0 0.0 0.0 ];

     elseif isequal(m,[1,2,1]) % m == [1,2,1]

         ordered_points(:,:,1,1) = [p1 p5 p8 p4 r4 r8 ze];
         ordered_points(:,:,1,4) = [p3 p7 p8 p4 r4 r8 ze];
         ordered_points(:,:,1,5) = [p2 p3 p4 p1 r4 r8 r5];
         ordered_points(:,:,1,6) = [p6 p7 p8 p5 r8 r4 r1];

         ordered_points(:,:,2,2) = [p2 p6 p7 p3 r3 r7 ze];
         ordered_points(:,:,2,4) = [p4 p8 p7 p3 r3 r7 ze];
         ordered_points(:,:,2,5) = [p1 p4 p3 p2 r3 r7 r6];
         ordered_points(:,:,2,6) = [p5 p8 p7 p6 r7 r3 r2];

         ordered_points(:,:,3,1) = [p1 p5 p8 p4 p3 p7 ze];
         ordered_points(:,:,3,2) = [p2 p6 p7 p3 p4 p8 ze];
         ordered_points(:,:,3,4) = [p4 p3 p7 p8 ze ze ze];
         ordered_points(:,:,3,5) = [p1 p2 p3 p4 p8 p7 ze];
         ordered_points(:,:,3,6) = [p5 p6 p7 p8 p4 p3 ze];

         ordered_points(:,:,5,1) = [p5 p8 p4 p1 r4 r3 r2];
         ordered_points(:,:,5,2) = [p6 p7 p3 p2 r3 r4 r1];
         ordered_points(:,:,5,4) = [p8 p7 p3 p4 r4 r3 ze];
         ordered_points(:,:,5,5) = [p1 p2 p3 p4 r4 r3 ze];

         ordered_points(:,:,6,1) = [p1 p4 p8 p5 r8 r7 r6];
         ordered_points(:,:,6,2) = [p2 p3 p7 p6 r7 r8 r5];
         ordered_points(:,:,6,4) = [p4 p3 p7 p8 r8 r7 ze];
         ordered_points(:,:,6,6) = [p5 p6 p7 p8 r8 r7 ze];


         adjacency_type = [ EAc 0.0 0.0 EAo VAo VAo
                            0.0 EAc 0.0 EAo VAo VAo
                            EAo EAo 0.0 ST  EAo EAo
                            0.0 0.0 0.0 0.0 0.0 0.0
                            VAo VAo 0.0 EAo EAc 0.0
                            VAo VAo 0.0 EAo 0.0 EAc ];

     elseif isequal(m,[2,2,1]) % m == [2,2,1]

         ordered_points(:,:,1,2) = [p2 p6 p7 p3 r4 r8 ze];
         ordered_points(:,:,1,4) = [p4 p8 p7 p3 r4 r8 ze];
         ordered_points(:,:,1,5) = [p1 p4 p3 p2 r4 r8 r5];
         ordered_points(:,:,1,6) = [p5 p8 p7 p6 r8 r4 r1];

         ordered_points(:,:,3,2) = [p2 p6 p7 p3 r2 r6 ze];
         ordered_points(:,:,3,4) = [p4 p8 p7 p3 r2 r6 ze];
         ordered_points(:,:,3,5) = [p1 p2 p3 p4 r2 r6 r5];
         ordered_points(:,:,3,6) = [p5 p6 p7 p8 r6 r2 r1];

         ordered_points(:,:,5,2) = [p6 p7 p3 p2 r4 r3 r2];
         ordered_points(:,:,5,4) = [p8 p7 p3 p4 r2 r3 r4];
         ordered_points(:,:,5,5) = [p1 p2 p3 p4 r2 r3 r4];

         ordered_points(:,:,6,2) = [p2 p3 p7 p6 r8 r7 r6];
         ordered_points(:,:,6,4) = [p4 p3 p7 p8 r6 r7 r8];
         ordered_points(:,:,6,6) = [p5 p6 p7 p8 r6 r7 r8];


         adjacency_type = [ 0.0 EAc 0.0 EAo VAo VAo
                            0.0 0.0 0.0 0.0 0.0 0.0
                            0.0 EAo 0.0 EAc VAo VAo
                            0.0 0.0 0.0 0.0 0.0 0.0
                            0.0 VAo 0.0 VAo VAc 0.0
                            0.0 VAo 0.0 VAo 0.0 VAc ];

     elseif isequal(m,[1,2,2]) % m == [1,2,2]

         ordered_points(:,:,1,1) = [p1 p4 p8 p5 r4 r8 r5];
         ordered_points(:,:,1,4) = [p3 p7 p8 p4 r5 r8 r4];
         ordered_points(:,:,1,6) = [p6 p7 p8 p5 r4 r8 r5];

         ordered_points(:,:,2,2) = [p2 p3 p7 p6 r3 r7 r6];
         ordered_points(:,:,2,4) = [p4 p8 p7 p3 r6 r7 r3];
         ordered_points(:,:,2,6) = [p5 p8 p7 p6 r3 r7 r6];

         ordered_points(:,:,3,1) = [p1 p5 p8 p4 r5 r6 r2];
         ordered_points(:,:,3,2) = [p2 p6 p7 p3 r6 r5 r1];
         ordered_points(:,:,3,4) = [p4 p3 p7 p8 r5 r6 ze];
         ordered_points(:,:,3,6) = [p5 p6 p7 p8 r5 r6 ze];

         ordered_points(:,:,5,1) = [p1 p4 p8 p5 r4 r3 r2];
         ordered_points(:,:,5,2) = [p2 p3 p7 p6 r3 r4 r1];
         ordered_points(:,:,5,4) = [p4 p3 p7 p8 r4 r3 ze];
         ordered_points(:,:,5,6) = [p5 p6 p7 p8 r4 r3 ze];


         adjacency_type = [ VAc 0.0 0.0 VAo 0.0 VAo
                            0.0 VAc 0.0 VAo 0.0 VAo
                            VAo VAo 0.0 EAc 0.0 EAo
                            0.0 0.0 0.0 0.0 0.0 0.0
                            VAo VAo 0.0 EAo 0.0 EAc
                            0.0 0.0 0.0 0.0 0.0 0.0 ]; 

     elseif isequal(m,[2,2,2]) % m == [2,2,2]

         ordered_points(:,:,1,2) = [p2 p3 p7 p6 r4 r8 r5];
         ordered_points(:,:,1,4) = [p4 p8 p7 p3 r5 r8 r4];
         ordered_points(:,:,1,6) = [p5 p8 p7 p6 r4 r8 r5];

         ordered_points(:,:,3,2) = [p2 p6 p7 p3 r5 r6 r2];
         ordered_points(:,:,3,4) = [p4 p3 p7 p8 r2 r6 r5];
         ordered_points(:,:,3,6) = [p5 p6 p7 p8 r2 r6 r5];

         ordered_points(:,:,5,2) = [p2 p3 p7 p6 r4 r3 r2];
         ordered_points(:,:,5,4) = [p4 p3 p7 p8 r2 r3 r4];
         ordered_points(:,:,5,6) = [p5 p6 p7 p8 r2 r3 r4];


         adjacency_type = [ 0.0 VAc 0.0 VAo 0.0 VAo
                            0.0 0.0 0.0 0.0 0.0 0.0
                            0.0 VAo 0.0 VAc 0.0 VAo
                            0.0 0.0 0.0 0.0 0.0 0.0
                            0.0 VAo 0.0 VAo 0.0 VAc
                            0.0 0.0 0.0 0.0 0.0 0.0 ];


    end

end



