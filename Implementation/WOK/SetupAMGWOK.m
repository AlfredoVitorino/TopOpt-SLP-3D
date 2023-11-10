function [An,P,M,R,perm,ngrids] = SetupAMGWOK(str,mg,kv,dens,p,emin,elem,mr)
% SetupAMGWOK performs the setup for algebraic multigrid when the global stiffness matrix is not explicitly generated. %
% INPUT: str - structure with the problem data. 
%        mg - structure that contains several parameters used by the multigrid method. 
%        kv - element stiffness matrix, converted into a column vector. 
%        dens - element densities vector (filtered).
%        p - penalty parameter of the SIMP model.
%        emin - Young's modulus of the void material.
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: An - system matrix of the coarsest grid.
%         P - cell containing the prolongation matrices. 
%         M - cell containing the smoother iteration matrices. 
%         R - Cholesky factor of the system matrix on the coarsest grid / or the preconditioner to use PCG on the coarsest grid. 
%         perm - approximate minimum degree permutation vector of the system matrix on the coarsest grid.
%         ngrids - number of grids. 
% ---------- %

ngrids = mg.ngrids; 
omega = mg.omega;
theta = mg.theta; 
optheta = mg.optheta;
nt = mg.nt; 
tolr = mg.tolr; 
tolq = mg.tolq;
itp = mg.itp; 
kappa = mg.kappa;
opchol = mg.opchol;
nmaxchol = mg.nmaxchol;
Anbatch = mg.Anbatch;

M = cell(ngrids-1,1); % Smoother matrices
P = cell(ngrids-1,1); % Prolongation matrices 

n = 1; 
while (n < ngrids)
    [M{n},nnode] = SmootherMatrix(P,str,kv,dens,p,emin,omega,n,Anbatch,elem,mr);
    
    % Test space %
    rng(0); 
    V = rand(nnode,nt);     
    V = TestSpaceWOK(M{n},nnode,nt,V,tolr,tolq,100000,omega,n,P,kv,dens,p,emin,str,elem,mr);
    % ---------- % 
    
    % Strong connections %
    [S,SC,lambda] = StrConnectionsWOK(V,nnode,theta,optheta,n,P,kv,dens,p,emin,str,elem,mr);
    % ---------- %

    % Coarsening %
    [C,F] = Coarsening(nnode,lambda,S,SC);
    % ---------- %

    % Prolongation %
    P{n} = ProlongationDPLS(nnode,S,V,C,F,itp,kappa);
    % ---------- %
    
    if((n == ngrids-1) && (opchol == 0))
        [~,sizeAn] = SmootherMatrix(P,str,kv,dens,p,emin,omega,n,Anbatch,elem,mr);
        if(sizeAn > nmaxchol)
            mg.ngrids = mg.ngrids + 1;  
            disp(['The size of the system for the coarsest grid is greater than nmaxchol. ngrids was increased to ',num2str(mg.ngrids)]);
        end
    end
    ngrids = mg.ngrids;
    
    n = n + 1; 
end

[An,~] = GlobalStiffMatrixWOK(P,str,kv,dens,p,emin,omega,ngrids,Anbatch,elem,mr);

perm = amd(An); % Approximate minimum degree permutation vector of A{ngrids}
if (mg.opchol == 1)
    if (mg.opcond == 1)        
        R = ichol(An(perm,perm), struct('diagcomp',0.1));
    else
        R = diag(diag(An));
    end
else
    R = chol(An(perm,perm));
end

end