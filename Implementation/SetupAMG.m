function [A,P,M,L,U,R,perm,ngrids] = SetupAMG(K,mg)
% SetupAMG performs the setup for algebraic multigrid. %
% INPUT: K - global stiffness matrix.
%        mg - structure that contains several parameters used by the multigrid method. 
% OUTPUT: A - cell containing the system matrices on each grid.
%         P - cell containing the prolongation matrices. 
%         M - cell containing the smoother iteration matrices. 
%         L - cell containing the lower triangular smoother iteration matrices (for SSOR).
%         U - cell containing the upper triangular smoother iteration matrices (for SSOR).
%         R - Cholesky factor of the system matrix on the coarsest grid / or the preconditioner to use PCG on the coarsest grid. 
%         perm - approximate minimum degree permutation vector of the system matrix on the coarsest grid.
%         ngrids - number of grids. 
% ---------- %

ngrids = mg.ngrids; 
smoother = mg.smoother;
omega = mg.omega;
theta = mg.theta; 
optheta = mg.optheta;
nt = mg.nt; 
tolr = mg.tolr; 
tolq = mg.tolq;
itp = mg.itp; 
kappa = mg.kappa;
opcond = mg.opcond; 
opchol = mg.opchol;
nmaxchol = mg.nmaxchol;

A = cell(ngrids,1);  % System matrices on each grid
A{1} = K;
P = cell(ngrids-1,1); % Prolongation matrices 
M = cell(ngrids-1,1); % Smoother iteration matrices 
L = cell(ngrids-1,1); % Lower triangular smoother iteration matrices (for SSOR).
U = cell(ngrids-1,1); % Upper triangular smoother iteration matrices (for SSOR).

n = 1; 
while (n < ngrids)
    nnode = size(A{n},1);
    
    % Smoother % 
    if(smoother == 0) % Jacobi
        M{n} = (1/omega)*diag(diag(A{n}));
    elseif(smoother == 1) % Gauss-Seidel/SOR
        M{n} = (1/omega)*diag(diag(A{n})) + tril(A{n},-1);
    elseif(smoother == 2) % SSOR 
        L{n} = (1/omega)*diag(diag(A{n})) + tril(A{n},-1); 
        U{n} = diag(diag(A{n}).^-1)*L{n}';
        L{n} = (omega/(2-omega))*L{n}; 
        M{n} = L{n}*U{n};
    end
    % ---------- %
    
    % Test space %
    rng(0); 
    V = rand(nnode,nt);     
    V = TestSpace(A{n},M{n},nnode,nt,V,tolr,tolq,100000,opcond);
    % ---------- % 
    
    % Strong connections %
    [S,SC,lambda] = StrConnections(A{n},V,nnode,theta,optheta);
    % ---------- %

    % Coarsening %
    [C,F] = Coarsening(nnode,lambda,S,SC);
    % ---------- %

    % Prolongation %
    P{n} = ProlongationDPLS(nnode,S,V,C,F,itp,kappa);
    % ---------- %
    
    A{n+1} = P{n}'*A{n}*P{n};
        
    if((n == ngrids-1) && (opchol == 0))
        if(size(A{n+1},1) > nmaxchol)
            mg.ngrids = mg.ngrids + 1;  
            disp(['The size of the system for the coarsest grid is greater than nmaxchol. ngrids was increased to ',num2str(mg.ngrids)]);
            L{n+1} = []; U{n+1} = []; 
        end
    end
    
    ngrids = mg.ngrids;
    n = n + 1; 
end

perm = amd(A{ngrids}); % Approximate minimum degree permutation vector of A{ngrids}
if (mg.opchol == 1)
    if (mg.opcond == 1)        
        R = ichol(A{ngrids}(perm,perm), struct('diagcomp',0.1));
    else
        R = diag(diag(A{ngrids}));
    end
else
    R = chol(A{ngrids}(perm,perm));
end

end