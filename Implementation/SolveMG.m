function [x,it] = SolveMG(x,f,Ac,P,M,L,U,R,perm,ngrids,cycle,smoother,smoothiter1,smoothiter2,tolMG,maxiterMG,opchol,opcond,pcgp)
% SolveMG applies a CG algorithm with multigrid preconditioner to solve the linear system. %
% INPUT: x - initial guess. 
%        f - right hand side vector. 
%        Ac - cell containing the system matrices on each grid.
%        P - cell containing the prolongation matrices. 
%        M - cell containing the smoother iteration matrices. 
%        L - cell containing the lower triangular smoother iteration matrices (for SSOR).
%        U - cell containing the upper triangular smoother iteration matrices (for SSOR).
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
%        pcpg - structure that contains several parameters used by the PCG algorithm.
% OUTPUT: x - approximated solution. 
%         it - number of iterations. 
% ---------- %

r = f - Ac{1}*x; 
nr = norm(r)/norm(f); % residual relative norm
it = 1; % number of iterations 
%fprintf('it: %4d | nr: %d\n', it, nr); 
while(nr > tolMG && it <= maxiterMG)
    if(cycle == 0)
        s = VCycle(Ac,r,zeros(length(r),1),smoothiter1,smoothiter2,ngrids,1,smoother,P,M,L,U,R,perm,opchol,opcond,pcgp); 
    elseif(cycle == 1)
        s = WCycle(Ac,r,zeros(length(r),1),smoothiter1,smoothiter2,ngrids,1,smoother,P,M,L,U,R,perm,opchol,opcond,pcgp); 
    elseif(cycle == 2)
        s = FMVCycle(Ac,r,smoothiter1,smoothiter2,ngrids,1,smoother,P,M,L,U,R,perm,opchol,opcond,pcgp); 
    end
    rho = r'*s;
    if(it == 1)
        d = s;
    else
        beta = rho/gamma2;
        d = beta*d + s;
    end
    q = Ac{1}*d;
    alpha = rho/(d'*q);
    x = x + alpha*d;
    r = r - alpha*q;
    gamma2 = rho;
    nr = norm(r)/norm(f);
    it = it + 1; 
    %fprintf('it: %4d | nr: %d\n', it, nr); 
end

end