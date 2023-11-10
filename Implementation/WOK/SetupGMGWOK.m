function [An,P,M,R,perm] = SetupGMGWOK(str,mg,kv,x,p,emin,elem,mr) 
% SetupGMGWOK performs the setup for geometric multigrid when the global stiffness matrix is not explicitly generated. %
% INPUT: str - structure with the problem data. 
%        mg - structure that contains several parameters used by the multigrid method.
%        kv - element stiffness matrix, converted into a column vector. 
%        x - element densities vector (filtered).
%        p - penalty parameter of the SIMP model.
%        emin - Young's modulus of the void material.
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: An - system matrix of the coarsest grid.
%         P - cell containing the prolongation matrices. 
%         M - cell containing the smoother iteration matrices. 
%         R - Cholesky factor of the system matrix on the coarsest grid / or the preconditioner to use PCG on the coarsest grid. 
%         perm - approximate minimum degree permutation vector of the system matrix of the coarsest grid.
% ---------- %

ngrids = mg.ngrids;
omega = mg.omega; 

% Grid parameters %
nelx = zeros(1,ngrids); % Number of elements in the length(x) direction 
nely = zeros(1,ngrids); % Number of elements in the height(y) direction
nelz = zeros(1,ngrids); % Number of elements in the width(z) direction
nelx(1) = str.nelx; nely(1) = str.nely; nelz(1) = str.nelz; 
for i = 2:ngrids
    nelx(i) = nelx(i-1)/2; 
    nely(i) = nely(i-1)/2; 
    nelz(i) = nelz(i-1)/2; 
end
nx = nelx*elem.deg+1; ny = nely*elem.deg+1; nz = nelz*elem.deg+1; % Number of nodes in each direction  
nnodes = (nx.*ny).*nz; % Total number of nodes on each grid

% Operators construction % 
freedofs = cell(ngrids,1); % Indexes of the free (without support) degrees of freedom on each grid
freedofs{1} = str.freedofs;
P = cell(ngrids-1,1);  % Prolongation matrices 
for i = 1:(ngrids-1)
    if(elem.type == 1 || elem.deg == 1) % Lagrange element
        [P{i},freedofs{i+1}] = Prolongation(nx(i),ny(i),nz(i),nx(i+1),ny(i+1),nz(i+1),nnodes(i+1),freedofs{i}); 
    elseif(elem.type == 2) % Serendipity element
        [P{i},freedofs{i+1}] = ProlongationS(nelx(i),nely(i),nelz(i),nelx(i+1),nely(i+1),nelz(i+1),freedofs{i},elem);
    end
end

[An,M] = GlobalStiffMatrixWOK(P,str,kv,x,p,emin,omega,ngrids,mg.Anbatch,elem,mr);

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
