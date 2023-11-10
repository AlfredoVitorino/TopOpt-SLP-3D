function [x,it] = SolveMGWOK(x,An,P,M,R,perm,ngrids,cycle,smoother,smoothiter1,smoothiter2,tolMG,maxiterMG,opchol,opcond,kv,dens,p,emin,pcgp,str,elem,mr)
% SolveMG applies a CG algorithm with multigrid preconditioner to solve the linear system, supposing that the global stiffness matrix K is not explicitly available. %
% INPUT: x - initial guess. 
%        An - system matrix of the coarsest grid.
%        P - cell containing the prolongation matrices. 
%        M - cell containing the smoother iteration matrices. 
%        R - Cholesky factor of the system matrix on the coarsest grid / or the preconditioner to use PCG on the coarsest grid. 
%        perm - approximate minimum degree permutation vector of the system matrix on the coarsest grid.
%        ngrids - number of grids. 
%        cycle - cycle type (0 - Vcycle, 1 - Wcycle, 2 - FullVCycle).
%        smoother - smoother (0 - Jacobi, 1 - Gauss-Seidel/SOR, 2 - SSOR).
%        smoothiter1 - number of pre-smooth iterations.
%        smoothiter2 - number of post-smooth iterations.
%        tolMG - tolerance for convergence. 
%        maxiterMG - maximum number of iterations.
%        opchol - option of what to do when the dimension of the system on the coarsest grid is greater than 'nmaxchol' (0 - increase ngrids, 1 - use PCG).
%        opcond - preconditioner option for PCG to solve the system on the coarsest grid (0 - diag, 1 - ichol).
%        kv - element stiffness matrix, converted into a column vector.
%        dens - element densities vector (filtered). 
%        p - penalty parameter for the SIMP model. 
%        emin - Young's modulus of the void material
%        pcpg - structure that contains several parameters used by the PCG algorithm.
%        str - structure with the problem data. 
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: x - approximated solution. 
%         it - number of iterations. 
% ---------- %

r = str.f - KelemProd(str,x,kv,dens,p,emin,elem,mr);
nr = norm(r)/norm(str.f); % residual relative norm
it = 1; % number of iterations 
%fprintf('it: %4d | nr: %d\n', it, nr); 
while(nr > tolMG && it <= maxiterMG)
    if(cycle == 0)
        s = VCycleWOK(An,r,zeros(length(r),1),smoothiter1,smoothiter2,ngrids,1,smoother,P,M,R,perm,opchol,opcond,kv,dens,p,emin,pcgp,str,elem,mr);
    elseif(cycle == 1)
        s = WCycleWOK(An,r,zeros(length(r),1),smoothiter1,smoothiter2,ngrids,1,smoother,P,M,R,perm,opchol,opcond,kv,dens,p,emin,pcgp,str,elem,mr);
    elseif(cycle == 2)
        s = FMVCycleWOK(An,r,smoothiter1,smoothiter2,ngrids,1,smoother,P,M,R,perm,opchol,opcond,kv,dens,p,emin,pcgp,str,elem,mr);
    end
    rho = r'*s;
    if(it == 1)
        d = s;
    else
        beta = rho/gamma2;
        d = beta*d + s;
    end
    q = KelemProd(str,d,kv,dens,p,emin,elem,mr);
    alpha = rho/(d'*q);
    x = x + alpha*d;
    r = r - alpha*q;
    gamma2 = rho;
    nr = norm(r)/norm(str.f);
    it = it + 1; 
    %fprintf('it: %4d | nr: %d\n', it, nr); 
end

end