% modified FlexGMRES with Righ Preconditioners, and standard restart
%
% Flex = can use a different preconditioner at each iteration,
% but only supports right preconditioning, see Y. Saad,
% "Iterative methods for sparse linear systems"
%
%   Version by J. Fernandez Villena
%              Computational Prototyping Group, RLE at MIT
%   Modification from GMRES to FlexGMRES by FastFieldSolvers S.r.L.
%
%   Can use double Righ Preconditioners
%   Can use function handle for A and preconditioners
%
% INPUT:  A         N-by-N matrix or function
%         b         rhs
%         restart   inner iterations before restart
%         tol       relative tolerance for the residue
%         maxit     maximum number of outer iterations (total is restart*maxit)
%         R1        first right preconditioner for A, so that A/R1 ~ EYE
%         R2        second right preconditioner for A, so that A/(R1*R2) ~ EYE
%         x0        initial guess
%
% OUTPUT: x         solution vector
%         flag      0 if converged, 1 if not
%         relres    final relative residue: norm(b - A*x)/norm(b) 
%         iter      vector with [current internal iterations, current external iterations]
%         resvec    vector containing norm of residual at each iteration of GMRES
%         
% EXAMPLE:  
%           A = sprand(1000,1000,0.05) +  spdiags(rand(1000,1),0,1000,1000);
%           b = rand(1000,1);
%           xtrue = A\b;
% 
%           [L,U] = luinc(A,0.01);
% 
%           % Approximation 1: no preconditioning
%           tic_gmres    = tic;
%           [x1,flag,relres,iter,resvec] = pgmres(A,b,30,1e-6,10);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 1: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
% 
%           % Approximation 3: right preconditioning
%           tic_gmres    = tic;
%           [x1,flag,relres,iter,resvec] = pgmres(A,b,30,1e-6,10,L,U);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 3: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
% 
%
%
% SUBFUNCTIONS: fpgmres_iter.m, iterchk, iterapp
%
%

function [x,flag,relres,iter,resvec] = fpgmres(A,b,restart,tol,maxit,R1,R2,x0)

% Initialize variables.
if(nargin < 2 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 3 || isempty(restart))
   restart = 50;
end
if(nargin < 4 || isempty(tol))
   tol = 1e-3;
end
if(nargin < 5 || isempty(maxit))
   maxit = 10;
end
if(nargin < 6 || isempty(R1))
   R1 = [];
end
if(nargin < 7 || isempty(R2))
   R2 = [];
end
if(nargin < 8 || isempty(x0))
   x0 = zeros(size(b));
end


% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);


% Calculate initial preconditioned residual.
% r = b - A*x0;
if (nnz(x0))
    r = b - iterapp('mtimes',afun,atype,afcnstr,x0);
else
    r = b;
end

% Calculate rhs norm
bnorm = norm(b);

% initialize the residue vector and other variables
resvec = zeros(restart*maxit,1);
it = 1;
outit = 0;
resvec(it) = norm(r)/bnorm;
x = x0;

% Restarted loop until convergence or maximum reached
while(resvec(it) > tol) && (outit < maxit)
  
    % increase external loop counter
    outit = outit+1;
    
   % Call the gmres to perform restart iterations
   [x,r,p,resvec_inner] = fpgmres_iter(A,x,r,restart,R1,R2,tol*bnorm);
   resvec(it+1:it+p) = resvec_inner/bnorm;
   [resvec(it+1:it+p)];
   if(it+p >= restart)
       disp(['   '])
       disp(['Achieved min norm of residual at each outer iter ::: ', num2str(min(resvec(it+1:it+p)))])
       disp(['   '])
   end
   it = it + p;
   
    % % % Calculate relative residual.
    % % res = b - iterapp('mtimes',afun,atype,afcnstr,x);
    % % fprintf(1,'restart %d, ||r|| = %g ( resvec %g, tol %g, reltol %g)\n' ,outit, norm(res), resvec(it), tol, tol*bnorm);
    
end


% residual, and iteration count
resvec = resvec(1:it);
iter = [p outit];

% check convergence
relres = norm(b - iterapp('mtimes',afun,atype,afcnstr,x))/norm(b);
if (relres <= tol ) % converged
    flag = 0;
else
    flag = 1;
end

disp(['   '])
% if(flag==0); disp('Gmres converged to the desired tolerance within MAXIT iterations.')
% elseif(flag==1);disp('Gmes iterated MAXIT times but did not converge.')
% %elseif(flag==2);disp('Preconditioner was ill-conditioned.')
% %elseif(flag==3);disp('Gmres stagnated (two consecutive iterates were the same).')
% end
  
%disp(['Achieved rel. residual NORM(B-A*X)/NORM(B) ::: ', num2str(relres)]);
disp(['# of total/inner/outer iterations  ::: ', num2str(length(resvec)),' / ',num2str(iter(1)),' / ',num2str(iter(2))]);
disp(['Achieved minimum residual norm NORM(B-A*X) ::: ',num2str(min(resvec))])



