function [Ftrue] = CorrectFunctionValue(strT,xstar,p,elem,emin,genK,opsolver,pcgp,mg,mr) 
% CorrectFunctionValue calculates the correct objective function value, solving the equilibrium linear system on the finest mesh. %
% INPUT: strT - structure with the problem data converted to the density grid with linear elements.
%        xstar - approximated optimal element densities vector.
%        p - penalty parameter for the SIMP model.
%        elem - structure with element characteristics (deg - polynomial degree, type - element type (1 = Lagrange, 2 = serendipity).
%        emin - Young's modulus of the void material.
%        genK - explicit generation of matrix K for multigrid solvers (true = yes, false = no).
%        opsolver - linear system solver option.
%        pcgp - structure that contains several parameters used by the pcg system solver.
%        mg - structure that contains several parameters used by the multigrid method.
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: Ftrue - corrected objective function value. 
% ---------- %

% Construct the stiffness matrix for the density mesh with linear elements of size 1 %
k = ElemStiffMatrix(strT);
kv = k(:);
if(genK)
    Dofs = DegreesOfFreedom(strT,elem);
    iK = kron(Dofs, ones(24,1))'; 
    iK = iK(:);
    jK = kron(Dofs, ones(1,24))'; 
    jK = jK(:);
    if(opsolver == 3 || opsolver == 4) % Using multigrid
        Ac = GlobalStiffMatrix(strT,kv,iK,jK,xstar,p,elem,emin,mr,[]); 
    else
        K = GlobalStiffMatrix(strT,kv,iK,jK,xstar,p,elem,emin,mr,[]); 
        perm = amd(K);
    end
end

% Solve the linear system %
u = zeros(3*strT.nnodes,1);
if(opsolver == 3) % Using AMG 
    if(genK)
        [Ac,P,M,L,U,R,perm,mg.ngrids] = SetupAMG(Ac,mg);
        [u(strT.freedofs),~] = SolveMG(u(strT.freedofs),strT.f,Ac,P,M,L,U,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,pcgp);
    else
        [An,P,M,R,perm,mg.ngrids] = SetupAMGWOK(strT,mg,kv,xstar,p,emin,elem,mr);
        [u(strT.freedofs),~] = SolveMGWOK(u(strT.freedofs),An,P,M,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,kv,xstar,p,emin,pcgp,strT,elem,mr);
    end
elseif(opsolver == 4) % Using GMG
    if(genK)
        [Ac,P,M,L,U,R,perm] = SetupGMG(Ac,strT,mg,elem);
        [u(strT.freedofs),~] = SolveMG(u(strT.freedofs),strT.f,Ac,P,M,L,U,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,pcgp);
    else
        [An,P,M,R,perm] = SetupGMGWOK(strT,mg,kv,xstar,p,emin,elem,mr);
        [u(strT.freedofs),~] = SolveMGWOK(u(strT.freedofs),An,P,M,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,kv,xstar,p,emin,pcgp,strT,elem,mr);
    end
elseif ((opsolver == 1) && (~genK))
    [u(strT.freedofs),~,~,~] = NodalDispWOK(strT,kv,xstar,p,emin,u(str.freedofs),pcgp,elem,mr);
else
    [u(strT.freedofs),~,~,~] = NodalDisp(K,strT.f,perm,opsolver,u(strT.freedofs),pcgp); 
end

Ftrue = strT.f'*u(strT.freedofs); % Correct function value 

end