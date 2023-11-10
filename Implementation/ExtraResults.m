function [exres] = ExtraResults(str,volfrac,p,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,xstarDsgn,xstar,exop)
% ExtraResults gathers additional results based on the solution obtained by the algorithm. %
% INPUT: str - structure with the problem data.
%        volfrac - maximum volume fraction of the domain that the structure can occupy.
%        p - penalty parameter for the SIMP model. 
%        emin - Young's modulus of the void material
%        opfilter - filter option (0 = no filter, 1 = weighted average density filter, 2 = average density filter).
%        rmin - filter radius.
%        volineq - type of the volume constraint (true = "less than or equal to", false = "equal to").
%        elem - structure with element characteristics (deg - polynomial degree, type - element type (1 = Lagrange, 2 = serendipity).
%        slp - structure that contains several parameters used by the SLP algorithm.
%        opsolver - linear system solver option.
%        pcgp - structure that contains several parameters used by the pcg system solver.
%        mg - structure that contains several parameters used by the multigrid method.
%        prj - structure that contains the parameters used for rounding the optimal densities to 0 or 1.
%        genK - explicit generation of the global stiffness matrix (true = yes, false = no).
%        mr - structure that contains parameters used by the multiresolution method.
%        strDens - structure with density mesh data for multiresolution. 
%        xstarDsgn - approximated optimal design variables vector.
%        xstar - approximated optimal element densities vector.
%        exop - options to calculate the extra results. 
% OUTPUT: exres - structure with the extra results.  
% ---------- %

% Converts the coarse mesh information, freedofs and loads vector to the density (fine) mesh with linear elements
if(mr.op) 
    strT = ConvertDofsMR(str,mr); 
else
    strT = str;
end
elem.deg = 1;
mr.op = false;
nelem = strT.nelem;
elem.fix = 0; 
% ---------- %

% Calculate the reference function value when multiresolution or an element degree greater than 1 is used
if(exop.correctF)  
    exres.Ftrue = CorrectFunctionValue(strT,xstar,p,elem,emin,genK,opsolver,pcgp,mg,mr);
end
% ---------- %

% Post optimize, solving the problem again on the density mesh, with the multiresolution solution as initial guess
if(exop.postopt)  
    [t.prefilter,t.filter,t.setupK,t.setupPrec,t.system,t.grad,t.LP,t.fix,t.proj,t.other,t.total] = deal(0.0);
    strT.fixedDens = strDens.fixedDens; strT.freeDens = strDens.freeDens; 
    rmin = rmin*(mr.n);
    k = ElemStiffMatrix(strT);
    kv = k(:);
    if(genK)
        Dofs = DegreesOfFreedom(strT,elem);
        [iK,jK,perm] = SetupStiffMatrix(strT,xstar,p,elem,Dofs,mr,strDens);
        [W,w] = ElemNeighborhood(strT,rmin,opfilter);
    else
        iK = []; jK = []; perm = []; Dofs = []; W = [];
        w = ElemNeighborhoodWOK(str,rmin,opfilter); 
    end
%     xstarDsgn = xstar;
%     xstar = ApplyFilter(xstarDsgn,W,w,opfilter,false);
    disp(' ');
    disp('Solving the problem again on the density grid, with the multiresolution solution as initial guess');
    [exres.xpost,~,~,~,~,~,exres.Fpost,~,~,~,~,~,~] = StrSLP(strT,volfrac,p,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,[],xstarDsgn,xstar,W,w,k,kv,iK,jK,perm,Dofs,t);
end
% ---------- %

% Round the solution, assigning 1 to the greatest densities until filling the volume and 0 to the others and calculate the function value
if(exop.roundsol)
    [~,I] = sort(xstar,'descend'); 
    vt = floor(volfrac*nelem); 
    xround = zeros(length(xstar),1); 
    xround(I(1:vt)) = 1.0; 
    exres.Fround = CorrectFunctionValue(strT,xround,p,elem,emin,genK,opsolver,pcgp,mg,mr); 
    exres.xround = xround;
end
% ---------- %

% Calculate the function value for the fully solid structure (with all densities equal to 1)
if(exop.solidsol)
    xsolid = ones(length(xstar),1); 
    exres.Fsolid = CorrectFunctionValue(strT,xsolid,p,elem,emin,genK,opsolver,pcgp,mg,mr);
end
% ---------- %

end