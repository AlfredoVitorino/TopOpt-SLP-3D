function [x] = Smooth(A,M,L,U,b,x0,niter,smoother)
% Smooth performs 'niter' iterations of the smoother on the linear system Ax = b, with initial guess x0. %
% INPUT: A - system matrix.
%        M - smoother iteration matrix. 
%        L - lower triangular smoother iteration matrix (for SSOR).
%        U - upper triangular smoother iteration matrix (for SSOR).
%        b - right hand side vector. 
%        x0 - initial guess. 
%        niter - number of iterations. 
%        smoother - smoother (0 - Jacobi, 1 - Gauss-Seidel/SOR, 2 - SSOR).
% OUTPUT: x - approximate solution.
% ---------- %

x = x0; 
if(smoother == 0 || smoother == 1) % Jacobi or Gauss-Seidel/SOR
    for i = 1:niter
        x = x + M\(b - A*x); 
    end
elseif(smoother == 2) % SSOR
    for i = 1:niter
        x = x + U\(L\(b - A*x)); 
    end
end
end