% Preconditioned FlexGMRES internal iteration routine
%
% Flex = can use a different preconditioner at each iteration,
% but only supports right preconditioning, see Y. Saad,
% "Iterative methods for sparse linear systems"
%
% Modification from GMRES to FlexGMRES by FastFieldSolvers S.r.L.
%
% Solves A z = r for z, then returns x + z
% Generates Arnoldi relation A V(:,1:m) = V(:,1:m+1) H
%
% INPUT:  A      N-by-N matrix
%         x      current solution vector
%         r      N-by-1 preconditioned residual vector
%         m      number of GMRES iterations to perform
%         R1     first right preconditioner for A
%         R2     second right preconditioner for A
%         tol    specifies the tolerance of the method
% OUTPUT: x      updated solution vector
%         r      preconditioned residual vector
%         k      number of GMRES iterations actually performed
%         resvec vector containing norm of residual at each iteration of GMRES
%
%

function [x,r,k,resvec] = fpgmres_iter(A,x,r,m,R1,R2,tol)

if(isempty(R1))
   existR1 = 0;
else
   existR1 = 1;
  [R1type,R1fun,R1fcnstr] = iterchk(R1);
end
if(isempty(R2))
   existR2 = 0;
else
   existR2 = 1;
   [R2type,R2fun,R2fcnstr] = iterchk(R2);
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);

% Preallocations
%
% Remark: if less than m iterations are required (overall)
% to reach convergence, we allocate much more memory than neede

% Preallocate and Initialize V
V = zeros(size(r,1),m+1);
V(:,1) = r / norm(r);
% Preallocate Z
Z = zeros(size(r,1),m);
% Preallocate H and resvec
H = zeros(m+1,m);
resvec = zeros(m,1);

for k = 1:m
    
    Z(:,k) = V(:,k);
    
   % Find Z(:,k) using right preconditioning if available.
   % R1*R2 approx A, then Z(:,k) = inv(R1*R2)*v = R2\(R1\v)
   % Z(:,k) = R2 \ ( R1 \ V(:,k) );
   % Note that using matrix type R1 and R2 does NOT require FlexGMRES,
   % as these are fixed. The expectation is that R1 and R2 (if present)
   % are function calls.
   if(existR1) % first right preconditioning
       if strcmp(R1type,'matrix')
           Z(:,k) = R1 \ Z(:,k);
       else
           Z(:,k) = iterapp('mtimes',R1fun,R1type,R1fcnstr,Z(:,k));
       end
   end
   if(existR2) % second right preconditioning
       if strcmp(R2type,'matrix')
           Z(:,k) = R2 \ Z(:,k);
       else
           Z(:,k) = iterapp('mtimes',R2fun,R2type,R2fcnstr,Z(:,k));
       end
   end

   
   % Apply the system operator
   % w = A*w;
   w = iterapp('mtimes',afun,atype,afcnstr,Z(:,k));
    
          
   % Create next column of V and H
   H(1:k,k) = V(:,1:k)' * w;
   w = w - V(:,1:k)*H(1:k,k);

   H(k+1,k) = norm(w);
   V(:,k+1) = w / H(k+1,k);

   % Initialize right hand side of least-squares system
   rhs = zeros(k+1,1);
   rhs(1) = norm(r);

   % Solve least squares system; Calculate residual norm
   Hc = H(1:k+1,1:k);
   y = Hc \ rhs;
   res = rhs - Hc * y;
   resvec(k) = norm(res);
   
   % check for early convergence
   if resvec(k) < tol
        break;
   end
   
end

% Calculate solution and residual
resvec = resvec(1:k);

% obtain update on solution
xupd = Z(:,1:k)*y;

% apply update x = x + Z(:,1:k) * y;
x = x + xupd;
r = V(:,1:k+1) * res;
