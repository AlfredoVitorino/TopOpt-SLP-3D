function [up,pcgit,tp,ts] = NodalDisp(K,f,perm,opsolver,uold,pcgp)
% NodalDisp calculates the nodal displacements vector. % 
% INPUT: K - global stiffness matrix.
%        f - nodal loads vector.
%        perm - approximate minimum degree permutation vector of K. 
%        opsolver - linear system solver option (0 = Cholesky factorization, 1 = CG w/ preconditioner being the diagonal of K, 2 = CG w/ preconditioner being the incomplete Cholesky factorization of K).
%        uold - previous vector u obtained on the last iteration, initial guess for pcg.
%        pcgp - structure that contains several parameters used by the pcg system solver.
% OUTPUT: up - nodal displacements vector (of the free nodes only).
%         pcgit - number of PCG iterations.
%         tp - time consumed to construct the preconditioner.
%         ts - time consumed to solve the linear system. 
% ---------- %

if(opsolver == 0) % Cholesky
    tic;
    R = chol(K(perm,perm));
    up = R\(R'\f(perm));
    up(perm) = up;
    ts = toc;
    pcgit = 0; tp = 0;
elseif(opsolver == 1) % Pcg - diagonal
    tic;
    B = diag(diag(K));
    tp = toc;
    tic;
    [up,~,~,pcgit] = pcg(K,f,pcgp.tolPCG,pcgp.maxiterPCG,B,[],uold);
    ts = toc;
elseif(opsolver == 2) % Pcg - incomplete Cholesky
    tic;
    B = ichol(K(perm,perm), struct('diagcomp',0.1));
    tp = toc;
    tic;
    [up,~,~,pcgit] = pcg(K(perm,perm),f(perm),pcgp.tolPCG,pcgp.maxiterPCG,B,B',uold(perm));
    up(perm) = up;
    ts = toc;
end

end