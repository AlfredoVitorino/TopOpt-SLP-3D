function [up,pcgit,tp,ts] = NodalDispWOK(str,kv,dens,p,emin,uold,pcgp,elem,mr)
% NodalDispWOK computes the nodal displacements vector without constructing the global stiffness matrix.
% INPUT: str - structure with the problem data. 
%        kv - element stiffness matrix, converted into a column vector.
%        dens - element densities vector (filtered). 
%        p - penalty parameter for the SIMP model. 
%        emin - Young's modulus of the void material
%        uold - previous vector u obtained on the last iteration, initial guess for pcg.
%        pcgp - structure that contains several parameters used by the pcg system solver.
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: up - nodal displacements vector (of the free nodes only).
%         pcgit - number of PCG iterations.
%         tp - time consumed to construct the preconditioner.
%         ts - time consumed to solve the linear system. 
% ---------- %

myKv = @(x) KelemProd(str,x,kv,dens,p,emin,elem,mr);
tic;
B = DiagPrec(str,kv,dens,p,emin,elem,mr);
tp = toc;
myMv = @(x) MvProd(x,B);
tic;
[up,~,~,pcgit] = pcg(myKv,str.f,pcgp.tolPCG,pcgp.maxiterPCG,myMv,[],uold);   
ts = toc;

end