function [x] = FMVCycleWOK(An,b,niter1,niter2,ngrids,nv,smoother,P,M,R,perm,opchol,opcond,kv,dens,p,emin,pcgp,str,elem,mr)
% FMVCycle applies the Full Multigrid V-Cycle recursively, without constructing the global stiffness matrix. %
% INPUT: An - system matrix of the coarsest grid. 
%        b - right hand side vector.
%        niter1 - number of pre-smooth iterations. 
%        niter2 - number of post-smooth iterations. 
%        ngrids - number of grids. 
%        nv - grid level. 
%        smoother - smoother (0 - Jacobi, 1 - Gauss-Seidel/SOR, 2 - SSOR). 
%        P - cell containing the prolongation matrices. 
%        M - cell containing the smoother iteration matrices. 
%        R - Cholesky factor of the system matrix on the coarsest grid / or the preconditioner to use PCG on the coarsest grid. 
%        perm - approximate minimum degree permutation vector of the system matrix on the coarsest grid.
%        opchol - option of what to do when the dimension of the system on the coarsest grid is greater than 'nmaxchol' (0 - increase ngrids, 1 - use PCG).
%        opcond - Preconditioner option for PCG to solve the system on the coarsest grid (0 - diag, 1 - ichol).
%        kv - element stiffness matrix, converted into a column vector.
%        dens - element densities vector (filtered). 
%        p - penalty parameter for the SIMP model. 
%        emin - Young's modulus of the void material
%        pcpg - structure that contains several parameters used by the PCG algorithm.
%        str - structure with the problem data.
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: x - approximated solution. 
% ---------- %

rH = P{nv}'*b; 
if(nv == ngrids-1)
    if (opchol == 1)
        if (opcond == 1)
            [x,~,~] = pcg(An(perm,perm),rH(perm),pcgp.tolPCG,pcgp.maxiterPCG,R,R');
            x(perm) = x;
        else
            [x,~,~] = pcg(An,rH,pcgp.tolPCG,pcgp.maxiterPCG,R);
        end
    else
        x = R\(R'\rH(perm)); 
        x(perm) = x;
    end
else
    x = FMVCycleWOK(An,rH,niter1,niter2,ngrids,nv+1,smoother,P,M,R,perm,opchol,opcond,kv,dens,p,emin,pcgp,str,elem,mr);
end
x = P{nv}*x; 
x = VCycleWOK(An,b,x,niter1,niter2,ngrids,nv,smoother,P,M,R,perm,opchol,opcond,kv,dens,p,emin,pcgp,str,elem,mr);

end