function [x] = VCycleWOK(An,b,x0,niter1,niter2,ngrids,nv,smoother,P,M,R,perm,opchol,opcond,kv,dens,p,emin,pcgp,str,elem,mr)
% VCycle applies the Multigrid V-Cycle recursively, without constructing the global stiffness matrix. %
% INPUT: An - system matrix of the coarsest grid. 
%        b - right hand side vector.
%        x0 - initial guess. 
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

x = SmoothWOK(M{nv},b,x0,niter1,smoother,nv,P,kv,dens,p,emin,str,elem,mr); 
rh = b - AciProd(x,str,nv,P,kv,dens,p,emin,elem,mr); %rh = b - A{nv}*x; 
rH = P{nv}'*rh;
if(nv == ngrids-1)
    if (opchol == 1)
        if (opcond == 1)
            [eH,~,~] = pcg(An(perm,perm),rH(perm),pcgp.tolPCG,pcgp.maxiterPCG,R,R');
            eH(perm) = eH;
        else
            [eH,~,~] = pcg(An,rH,pcgp.tolPCG,pcgp.maxiterPCG,R);
        end
    else
        eH = R\(R'\rH(perm)); 
        eH(perm) = eH;
    end
else
    eH = VCycleWOK(An,rH,zeros(length(rH),1),niter1,niter2,ngrids,nv+1,smoother,P,M,R,perm,opchol,opcond,kv,dens,p,emin,pcgp,str,elem,mr); 
end
eh = P{nv}*eH;
x = x + eh; 
x = SmoothWOK(M{nv},b,x,niter2,smoother,nv,P,kv,dens,p,emin,str,elem,mr);

end