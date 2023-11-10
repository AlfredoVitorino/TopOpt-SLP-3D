function [x] = FMVCycle(A,b,niter1,niter2,ngrids,nv,smoother,P,M,L,U,R,perm,opchol,opcond,pcgp)
% FMVCycle applies the Full Multigrid V-Cycle recursively. %
% INPUT: A - cell containing the system matrices on each grid. 
%        b - right hand side vector.
%        niter1 - number of pre-smooth iterations. 
%        niter2 - number of post-smooth iterations. 
%        ngrids - number of grids. 
%        nv - grid level. 
%        smoother - smoother (0 - Jacobi, 1 - Gauss-Seidel/SOR, 2 - SSOR). 
%        P - cell containing the prolongation matrices. 
%        M - cell containing the smoother iteration matrices. 
%        L - cell containing the lower triangular smoother iteration matrices (for SSOR).
%        U - cell containing the upper triangular smoother iteration matrices (for SSOR).
%        R - Cholesky factor of the system matrix on the coarsest grid / or the preconditioner to use PCG on the coarsest grid. 
%        perm - approximate minimum degree permutation vector of the system matrix on the coarsest grid.
%        opchol - option of what to do when the dimension of the system on the coarsest grid is greater than 'nmaxchol' (0 - increase ngrids, 1 - use PCG).
%        opcond - Preconditioner option for PCG to solve the system on the coarsest grid (0 - diag, 1 - ichol).
%        pcpg - structure that contains several parameters used by the PCG algorithm.
% OUTPUT: x - approximated solution. 
% ---------- %

rH = P{nv}'*b; 
if(nv == ngrids-1)
    if (opchol == 1)
        if (opcond == 1)
            [x,~,~] = pcg(A{nv+1}(perm,perm),rH(perm),pcgp.tolPCG,pcgp.maxiterPCG,R,R');
            x(perm) = x;
        else
            [x,~,~] = pcg(A{nv+1},rH,pcgp.tolPCG,pcgp.maxiterPCG,R);
        end
    else
        x = R\(R'\rH(perm)); 
        x(perm) = x;
    end
else
    x = FMVCycle(A,rH,niter1,niter2,ngrids,nv+1,smoother,P,M,L,U,R,perm,opchol,opcond,pcgp);
end
x = P{nv}*x; 
x = VCycle(A,b,x,niter1,niter2,ngrids,nv,smoother,P,M,L,U,R,perm,opchol,opcond,pcgp);

end