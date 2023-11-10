function [x] = SmoothWOK(M,b,x0,niter,smoother,nv,P,kv,dens,p,emin,str,elem,mr)
% Smooth performs 'niter' iterations of the smoother on the linear system Ax = b, with initial guess x0, when the system matrix A is not explicitly available. %
% INPUT: M - smoother iteration matrix. 
%        b - right hand side vector. 
%        x0 - initial guess. 
%        niter - number of iterations. 
%        smoother - smoother (0 - Jacobi, 1 - Gauss-Seidel/SOR, 2 - SSOR).
%        nv - grid level. 
%        P - cell containing the prolongation matrices. 
%        kv - element stiffness matrix, converted into a column vector.
%        dens - element densities vector (filtered). 
%        p - penalty parameter for the SIMP model. 
%        emin - Young's modulus of the void material.
%        str - structure with the problem data. 
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: x - approximate solution.
% ---------- %

x = x0; 
if(smoother == 0 || smoother == 1) % Jacobi or Gauss-Seidel/SOR
    for i = 1:niter
        x = x + M\(b - AciProd(x,str,nv,P,kv,dens,p,emin,elem,mr));
    end
elseif(smoother == 2) % SSOR
    error('Error: mg.smoother must be 0 if genK = false.');
end

end