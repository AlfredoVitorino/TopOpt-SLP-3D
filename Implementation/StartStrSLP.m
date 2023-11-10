function [xstar,xstarDsgn,u,opstop,iter,itrej,F,vol,ns,pcgit,lpit,time,xhist,strDens] = StartStrSLP(str,volfrac,p,p123,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr)
% StartStrSLP initializes the SLP algorithm to solve the topology optimization problem. %
% INPUT: str - structure with the problem data.
%        volfrac - maximum volume fraction of the domain that the structure can occupy.
%        p - penalty parameter for the SIMP model. 
%        p123 - indicates if the continuation strategy (use p = 1,2,3) will be adopted. 
%        emin - Young's modulus of the void material.
%        opfilter - filter option (0 = no filter, 1 = weighted average density filter, 2 = average density filter).
%        rmin - filter radius.
%        volineq - type of the volume constraint (true = "less than or equal to", false = "equal to").
%        elem - structure with element characteristics (deg - polynomial degree, type - element type (1 = Lagrange, 2 = serendipity).
%        slp - structure that contains several parameters used by the SLP algorithm.
%        opsolver - linear system solver option.
%        pcgp - structure that contains several parameters used by the pcg system solver.
%        mg - structure that contains several parameters used by the multigrid method.
%        prj - structure that contains the parameters used for rounding the optimal densities to 0 or 1.
%        genK - explicit generation of matrix K for multigrid solvers (true = yes, false = no).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: xstar - approximated optimal element densities vector. 
%         xstarDsgn - approximated optimal design variables vector.
%         u - approximated nodal displacements vector. 
%         opstop - stopping criterion achieved.
%         iter - number of iterations.
%         itrej - number of rejected steps.
%         F - objective function optimal value.
%         vol - volume fraction at the optimal solution.
%         ns - number of linear systems solved.
%         pcgit - number of iterations of PCG for each system solved.
%         lpit - number of iterations of linprog for each LP solved.
%         time - structure with the time consumed on each part of the program.
%         xhist - densities of each iteration of the SLP algorithm.
%         strDens - structure with density mesh data for multiresolution. 
% ---------- %

disp('Solving the topology optimization problem using a Sequential Linear Programming algorithm.')

startTime = tic;

t.prefilter = 0.0; % Time consumed to obtain the element neighbors and weight factors
t.filter = 0.0; % Time consumed to apply the filter
t.setupK = 0.0; % Time consumed to construct the stiffness matrices 
t.setupPrec = 0.0; % Time consumed to setup the preconditioner to solve the linear systems
t.system = 0.0; % Time consumed to solve the linear systems
t.grad = 0.0; % Time consumed to calculate the gradient of the objective function
t.LP = 0.0; % Time consumed to solve the LP subproblems 
t.fix = 0.0; % Time consumed to choose and fix elements
t.proj = 0.0; % Time consumed to round the densities to 0 or 1 

% Checking if there are void or solid regions in the domain %
if(~isfield(str,'fixedDens'))
    str.freeDens = 1:str.nelem;
    str.fixedDens = []; 
    str.fixedDensVal = []; 
end
% ---------- %

% Getting the initial design variables and meshes data for multiresolution %
if(~mr.op)   
    strDens = [];
    strDsgn = [];
    x0 = volfrac*ones(str.nelem,1);
    if(~isempty(str.fixedDens)) % Assign fixed design variables
        x0(str.fixedDens) = str.fixedDensVal;
        x0(str.freeDens) = (volfrac*str.l*str.h*str.w-sum(str.fixedDensVal)*str.el*str.eh*str.ew)/length(str.freeDens);
    end
else 
    [strDens,strDsgn,x0] = StartMR(str,volfrac,p,p123,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr);
end
% ---------- %

% Calculating the matrix to interpolate the displacements %
if(mr.op && mr.interp)
    [mr.Idisp, mr.Ldofs, mr.kd] = InterpolateDisp(mr,elem,strDens);
end
% ---------- %

% Adapting the number of nodes, freedofs, vector of nodal loads and multigrid parameters according to the polynomial degree of the shape functions %
if(elem.deg > 1)
    str = ConvertDofs(str,elem);
    if(genK)
        mg.smoother = 2; mg.omega = 1.2;
    else
        mg.smoother = 0; mg.omega = 0.3;
    end
end
% ---------- %

% Disabling the incomplete Cholesky preconditioner and warning that genK = false can only be used with opsolver = 1, 3 or 4 %
if (~genK)
    if (opsolver == 4 || opsolver == 3)
        if (mg.smoother ~= 0)
            disp('Using Multigrid with genK = false. Smoother will be set to 0 (Jacobi).');
            mg.smoother = 0;
        end
        if(opsolver == 3 && mg.opcond ~= 0) 
            disp('Using Algebraic Multigrid with genK = false. The preconditioner for test space generation will be set to 0 (diagonal).');
            mg.opcond = 1; 
        end
    elseif (opsolver ~= 1)
        disp(['opsolver = ',num2str(opsolver),'. genK will be set to true.']);
        genK = true;
    end
end
% ---------- %

% Checking if the geometric multigrid can be used (if not, try to reduce ngrids or change to algebraic multigrid) %
if (opsolver == 4 && (mod(str.nelx,2^(mg.ngrids-1))+mod(str.nely,2^(mg.ngrids-1))+mod(str.nelz,2^(mg.ngrids-1))) ~= 0)
    ngridstemp = mg.ngrids;
    while((mod(str.nelx,2^(mg.ngrids-1))+mod(str.nely,2^(mg.ngrids-1))+mod(str.nelz,2^(mg.ngrids-1))) ~= 0)
        mg.ngrids = mg.ngrids - 1;
    end
    if (mg.ngrids >= 2)
        disp(['You cannot use geometric multigrid in this problem with ',num2str(ngridstemp),' grids. ngrids was reduced to ',num2str(mg.ngrids)]);
    else
        disp('You cannot use geometric multigrid in this problem. Using algebraic multigrid instead.');
        opsolver = 3;
        mg.ngrids = ngridstemp;
    end        
end
% ---------- %

% Adjusting ngrids and opchol according to nmaxchol (for geometric multigrid) %
if ((opsolver == 4) && (mg.opchol == 0)) 
    while ((((str.nelx/(2^(mg.ngrids-1))+1)*(str.nely/(2^(mg.ngrids-1))+1)*(str.nelz/(2^(mg.ngrids-1))+1)*3)>mg.nmaxchol)&&(mg.opchol == 0))
        if ((mod(str.nelx,2^mg.ngrids)+mod(str.nely,2^mg.ngrids)+mod(str.nelz,2^mg.ngrids)) == 0)
            mg.ngrids = mg.ngrids+1;
            disp(['The size of the system for the coarsest grid is greater than nmaxchol. ngrids was increased to ',num2str(mg.ngrids)]);
        else
            mg.opchol = 1;
            disp('The size of the system for the coarsest grid is greater than nmaxchol. PCG will be used.');
        end
    end
end
% ---------- %

% Element neighbors and weight factors (for the filter application) %
if(opfilter ~= 0)
    tic;
    if(genK)
        if(mr.op) % Multiresolution          
            if(mr.n ~= mr.d) % Multiresolution with 3 distinct meshes
                [W,w] = ElemNeighborhoodMR(strDens,strDsgn,rmin,mr,opfilter);
            else
                [W,w] = ElemNeighborhood(strDsgn,rmin,opfilter);
            end
        else
            [W,w] = ElemNeighborhood(str,rmin,opfilter);
        end
    else
        if(mr.op)
            if(mr.n ~= mr.d)
                w = ElemNeighborhoodMRWOK(strDens,strDsgn,rmin,mr,opfilter); 
            else
                w = ElemNeighborhoodWOK(strDsgn,rmin,opfilter); 
            end
        else
            w = ElemNeighborhoodWOK(str,rmin,opfilter); 
        end
        W = [];
    end
    t.prefilter = t.prefilter + toc; 
else
    W = []; w = [];
end
% ---------- % 

% Setting up finite element analysis %
tic;
if(mr.op) % Multiresolution
    [k,kv] = ElemStiffMatrixGQMR(str,mr.n,elem); % Stiffness integrands for each density element inside a finite element
else
    if(elem.deg > 1)
        k = ElemStiffMatrixGQ(str,elem.deg+1,elem); % Element stiffness matrix with Gaussian quadrature
    else
        k = ElemStiffMatrix(str); % Element stiffness matrix
    end
    kv = k(:); % Element stiffness matrix converted into a vector
end
if(genK)
    Dofs = DegreesOfFreedom(str,elem); 
    xfil0 = ApplyFilter(x0,W,w,opfilter,false);
    [iK,jK,perm] = SetupStiffMatrix(str,xfil0,p,elem,Dofs,mr,strDens);
else
    iK = []; jK = []; perm = []; Dofs = []; 
    xfil0 = ApplyFilterWOK(x0,str,rmin,w,opfilter,false,mr,strDens,strDsgn);
end
t.setupK = t.setupK + toc; 
% ---------- %

% Disabling some parameters in the first call to StrSLP %
if(~prj.op)
    temp.heavi = false; % Disabling the Heaviside projection if the densities are not to be rounded to 0 or 1.
else
    temp.heavi = prj.heavi; % Storing prj.heavi.
end
prj.heavi = false; % Disabling the Heaviside projection in the first call to StrSLP.
temp.fix = elem.fix; % Storing elem.fix. 
elem.fix = 0; % Don't fix variables in the first call to StrSLP.
temp.prjop = prj.op; % Storing prj.op
temp.tolG = slp.tolG; % Storing slp.tolG
if(elem.increasedeg) 
    prj.op = false; % Don't round densities in the first call to StrSLP if we are increasing the element degree.  
    slp.tolG = 10*slp.tolG; % Increase the minimum value for the projected gradient norm in the first call to StrSLP if we are increasing the element degree. 
end
% ---------- %

% Applying the SLP algorithm to solve the problem %
if (~p123)  
    [xstar,xstarDsgn,u,opstop,iter,itrej,F,vol,ns,pcgit,lpit,t,xhist] = StrSLP(str,volfrac,p,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,strDsgn,x0,xfil0,W,w,k,kv,iK,jK,perm,Dofs,t); 
else
    tmp = slp.maxiter;  
    slp.maxiter = slp.maxit12;
    tprj = prj.op;
    prj.op = false;
    disp('Solving the problem with SIMP parameter p = 1.');
    [xstar,xstarDsgn,~,~,iter1,itrej1,~,~,ns1,pcgit1,lpit1,t,~] = StrSLP(str,volfrac,1,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,strDsgn,x0,xfil0,W,w,k,kv,iK,jK,perm,Dofs,t); 
    disp('Solving the problem with SIMP parameter p = 2.');
    [xstar,xstarDsgn,~,~,iter2,itrej2,~,~,ns2,pcgit2,lpit2,t,~] = StrSLP(str,volfrac,2,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,strDsgn,xstarDsgn,xstar,W,w,k,kv,iK,jK,perm,Dofs,t); 
    slp.maxiter = tmp;
    prj.op = tprj;
    disp('Solving the problem with SIMP parameter p = 3.');
    [xstar,xstarDsgn,u,opstop,iter3,itrej3,F,vol,ns3,pcgit3,lpit3,t,xhist] = StrSLP(str,volfrac,3,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,strDsgn,xstarDsgn,xstar,W,w,k,kv,iK,jK,perm,Dofs,t);  
    iter = [iter1,iter2,iter3];
    itrej = [itrej1,itrej2,itrej3];
    ns = [ns1,ns2,ns3];
    pcgit = [pcgit1;pcgit2;pcgit3];
    lpit = [lpit1;lpit2;lpit3];
end
% ---------- %

% Calculating the time consumed %
total = toc(startTime); % Total time
t.other = total - (t.prefilter + t.filter + t.setupK + t.setupPrec + t.system + t.grad + t.LP + t.fix + t.proj);
t.total = total;
temp.time = struct2table(t); 
temp.rownames = "First solution";
% ---------- %
 
% Storing the number of iterations %
temp.iter = iter;
temp.itrej = itrej;
temp.ns = ns;
temp.pcgit = pcgit;
temp.lpit = lpit;
% ---------- %

% Solve the problem again if we are increasing the element degree %
if (elem.increasedeg)   
    if (genK) % Change smoother parameters to improve multigrid with element degree greater than 1
        mg.smoother = 2; mg.omega = 1.2;
    else
        mg.smoother = 0; mg.omega = 0.3;
    end
    while (elem.deg < elem.maxdeg) 
        % Restart the timer
        startTime = tic;
        [t.prefilter,t.filter,t.setupK,t.setupPrec,t.system,t.grad,t.LP,t.fix,t.proj,t.other,t.total] = deal(0.0);

        elem.fix = temp.fix; 
        x0 = xstarDsgn; 
        xfil0 = xstar;
        elem.deg = elem.deg + 1;
        oldstr = str;

        % Use serendipity elements if the degree is greater than 2
        if(elem.deg > 2)
            elem.type = 2;
        end

        % Recover the defined minimum value for the projected gradient norm if the element degree is maximum 
        if(elem.deg == elem.maxdeg)
            slp.tolG = temp.tolG; 
            prj.op = temp.prjop;
        end

        % Adapting the number of nodes, freedofs and vector of nodal loads according to the polynomial degree of the shape functions
        str = ConvertDofs(str,elem);

        % Saving the original freedofs and vector of nodal loads if we will choose variables to be fixed 
        if(elem.fix ~= 0)
            str.oldfreedofs = str.freedofs;
            str.suppindex = setdiff(1:3*str.nnodes,str.oldfreedofs);
            str.oldf = zeros(3*str.nnodes,1);
            str.oldf(str.oldfreedofs) = str.f;          
        end

        % Recalculating the matrix to interpolate the displacements 
        if(mr.op && mr.interp)
            [mr.Idisp, mr.Ldofs, mr.kd] = InterpolateDisp(mr,elem,strDens);
        end

        % Setting up finite element analysis 
        tic;
        if(mr.op) % Multiresolution
            [k,kv] = ElemStiffMatrixGQMR(str,mr.n,elem); % Stiffness integrands for each density element inside a finite element
        else
            if(elem.deg > 1)
                k = ElemStiffMatrixGQ(str,elem.deg+1,elem); % Element stiffness matrix with Gaussian quadrature
            else
                k = ElemStiffMatrix(str); % Element stiffness matrix
            end
            kv = k(:); % Element stiffness matrix converted into a vector
        end
        if(genK)
            Dofs = DegreesOfFreedom(str,elem); 
            [iK,jK,perm] = SetupStiffMatrix(str,xstar,p,elem,Dofs,mr,strDens);
        end
        t.setupK = t.setupK + toc; 

        % Calling StrSLP
        disp(' ');
        disp(['Solving the problem with element degree equal to ', num2str(elem.deg)]);
        [xstar,xstarDsgn,u,opstop,iter4,itrej4,F,vol,ns4,pcgit4,lpit4,t,xhist] = StrSLP(str,volfrac,p,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,strDsgn,x0,xfil0,W,w,k,kv,iK,jK,perm,Dofs,t);

        % Storing the number of iterations 
        temp.iter = [temp.iter, iter4];
        temp.itrej = [temp.itrej, itrej4];
        temp.ns = [temp.ns, ns4];
        temp.pcgit = [temp.pcgit; pcgit4];
        temp.lpit = [temp.lpit; lpit4];

        % Calculating the time consumed
        total = toc(startTime); % Total time
        t.other = total - (t.prefilter + t.filter + t.setupK + t.setupPrec + t.system + t.grad + t.LP + t.fix + t.proj);
        t.total = total;
        temp.time = [temp.time; struct2table(t)];
        temp.rownames = [temp.rownames; "Solution w/ deg "+num2str(elem.deg)]; 
    end
end
% ---------- %

% Solve the problem again when rounding the densities, until find an optimal solution %
prj.op = temp.prjop;
if (prj.op)
    % Restart the timer
    startTime = tic;
    [t.prefilter,t.filter,t.setupK,t.setupPrec,t.system,t.grad,t.LP,t.fix,t.proj,t.other,t.total] = deal(0.0);

    % Restoring some values 
    prj.heavi = temp.heavi; 
    opfilter = prj.opfilter;

    % Saving the original freedofs and vector of nodal loads if we will choose elements to be fixed 
    if(elem.fix ~= 0)
        str.oldfreedofs = str.freedofs;
        str.suppindex = setdiff(1:3*str.nnodes,str.oldfreedofs);
        str.oldf = zeros(3*str.nnodes,1);
        str.oldf(str.oldfreedofs) = str.f;
    end
    
    % Use linear elements if prj.deg = true
    if(elem.deg > 1 && prj.deg)
        elem.deg = 1; 
        str = oldstr;

        % Recalculating the matrix to interpolate the displacements 
        if(mr.op && mr.interp)
            [mr.Idisp, mr.Ldofs] = InterpolateDisp(mr,elem);
        end

        % Setting up finite element analysis 
        tic;
        if(mr.op) % Multiresolution
            [k,kv] = ElemStiffMatrixGQMR(str,mr.n,elem); % Stiffness integrands for each density element inside a finite element
        else
            if(elem.deg > 1)
                k = ElemStiffMatrixGQ(str,elem.deg+1,elem); % Element stiffness matrix with Gaussian quadrature
            else
                k = ElemStiffMatrix(str); % Element stiffness matrix
            end
            kv = k(:); % Element stiffness matrix converted into a vector
        end
        if (genK)
            Dofs = DegreesOfFreedom(str,elem); 
            [iK,jK,perm] = SetupStiffMatrix(str,xstar,p,elem,Dofs,mr,strDens);
        end
        t.setupK = t.setupK + toc; 
    end

    % Recomputing element neighbors and weight factors (for the filter application)      
    if(opfilter ~= 0)
        if (prj.rmin < rmin && ~elem.increasedeg)
            rmin = prj.rmin;
            tic;
            if(genK)
                if(mr.op) % Multiresolution          
                    if(mr.n ~= mr.d) % Multiresolution with 3 distinct meshes
                        [W,w] = ElemNeighborhoodMR(strDens,strDsgn,rmin,mr,opfilter);
                    else
                        [W,w] = ElemNeighborhood(strDsgn,rmin,opfilter);
                    end
                else
                    [W,w] = ElemNeighborhood(str,rmin,opfilter);
                end
            else
                if(mr.op)
                    if(mr.n ~= mr.d)
                        w = ElemNeighborhoodMRWOK(strDens,strDsgn,rmin,mr,opfilter); 
                    else
                        w = ElemNeighborhoodWOK(strDsgn,rmin,opfilter); 
                    end
                else
                    w = ElemNeighborhoodWOK(str,rmin,opfilter); 
                end
                W = [];
            end
            t.prefilter = t.prefilter + toc; 
        end
    else
        W = []; w = [];
    end

    tic;
    if (prj.heavi)
        if(mr.op)
            velem = ((str.l*str.h*str.w)/strDens.nelem)*ones(strDens.nelem,1);
        else
            velem = ((str.l*str.h*str.w)/str.nelem)*ones(str.nelem,1);
        end 
        vcur = velem'*xstar;
        vdes = str.l*str.h*str.w*volfrac;
        [xstar, prj] = HeavisideProj(xstar,prj,velem,vcur,vdes,true);
    end

    rndfrac = round(length(xstar)*volfrac)/length(xstar);
    t.proj = t.proj + toc;
    eleq0 = length(find(xstar<=0));
    eleq1 = length(find(xstar>=1));
    elint = length(xstar)-eleq0-eleq1;
    disp(['Void elements: ',num2str(eleq0),'.  Full elements: ',num2str(eleq1),'. Intermediate elements: ',num2str(elint),'. Volume fraction: ',num2str(vol*100)]);

    % Calling StrSLP %
    itproj = 0;
    xold = zeros(length(xstar),1);
    while ((itproj<prj.maxit)&&((norm(xold-xstar,1)>(norm(xold,1)*prj.tolN))||((vol-rndfrac)>prj.tolV)))
        xold = xstar;
        if ((mr.op)&&(mr.n~=mr.d))
            x0 = xstarDsgn;
        else
            x0 = xold;
        end       
        
        tic;
        if (genK)
            xfil0 = ApplyFilter(x0,W,w,opfilter,false);
        else
            xfil0 = ApplyFilterWOK(x0,str,rmin,w,opfilter,false,mr,strDens,strDsgn);
        end
        t.filter = t.filter + toc;
        
        disp(' ');
        disp(['Attempt ',num2str(itproj+1),' to solve the problem after rounding the densities to 0 or 1 using the gradient of the Lagrangian.'])
        [xstar,xstarDsgn,uold,opstop,iter4,itrej4,~,vol,ns4,pcgit4,lpit4,t,~] = StrSLP(str,volfrac,p,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,strDsgn,x0,xfil0,W,w,k,kv,iK,jK,perm,Dofs,t); 
        
        if (prj.heavi)
            tic;
            [xstar, prj] = HeavisideProj(xstar,prj,velem,velem'*xstar,vdes,true);
            t.proj = t.proj + toc;
        end
        eleq0 = length(find(xstar<=0));
        eleq1 = length(find(xstar>=1));
        elint = length(xstar)-eleq0-eleq1;
        disp(['Void elements: ',num2str(eleq0),'.  Full elements: ',num2str(eleq1),'. Intermediate elements: ',num2str(elint),'. Volume fraction: ',num2str(vol*100)]);
        itproj = itproj+1;
        
        % Storing the number of iterations 
        temp.iter = [temp.iter, iter4];
        temp.itrej = [temp.itrej, itrej4];
        temp.ns = [temp.ns, ns4];
        temp.pcgit = [temp.pcgit; pcgit4];
        temp.lpit = [temp.lpit; lpit4]; 
    end

    % Calculating the time consumed
    total = toc(startTime); % Total time
    t.other = total - (t.prefilter + t.filter + t.setupK + t.setupPrec + t.system + t.grad + t.LP + t.fix + t.proj);
    t.total = total;
    temp.time = [temp.time; struct2table(t)]; 
    temp.rownames = [temp.rownames; "Round densities"]; 

    % Calculate the objective function value of the final solution with rounded densities 
    u = zeros(3*str.nnodes,1); 
    if(opsolver == 3) % Using AMG 
        if(genK)
            Ac = GlobalStiffMatrix(str,kv,iK,jK,xstar,p,elem,emin,mr,strDens);
            [Ac,P,M,L,U,R,perm,mg.ngrids] = SetupAMG(Ac,mg);
            [u(str.freedofs),~] = SolveMG(uold(str.freedofs),str.f,Ac,P,M,L,U,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,pcgp);
        else
            [An,P,M,R,perm,mg.ngrids] = SetupAMGWOK(str,mg,kv,xstar,p,emin,elem,mr);
            [u(str.freedofs),~] = SolveMGWOK(uold(str.freedofs),An,P,M,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,kv,xstar,p,emin,pcgp,str,elem,mr);
        end
    elseif(opsolver == 4) % Using GMG
        if(genK)
            Ac = GlobalStiffMatrix(str,kv,iK,jK,xstar,p,elem,emin,mr,strDens);
            [Ac,P,M,L,U,R,perm] = SetupGMG(Ac,str,mg,elem);
            [u(str.freedofs),~] = SolveMG(uold(str.freedofs),str.f,Ac,P,M,L,U,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,pcgp);
        else
            [An,P,M,R,perm] = SetupGMGWOK(str,mg,kv,xstar,p,emin,elem,mr);
            [u(str.freedofs),~] = SolveMGWOK(uold(str.freedofs),An,P,M,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,kv,xstar,p,emin,pcgp,str,elem,mr);
        end
    elseif ((opsolver == 1) && (~genK))
        [u(str.freedofs),~,~,~] = NodalDispWOK(str,kv,xstar,p,emin,uold(str.freedofs),pcgp,elem,mr);
    else
        K = GlobalStiffMatrix(str,kv,iK,jK,xstar,p,elem,emin,mr,strDens); 
        [u(str.freedofs),~,~,~] = NodalDisp(K,str.f,perm,opsolver,uold(str.freedofs),pcgp); 
    end
    F = str.f'*u(str.freedofs);
end
% ---------- %

% Total number of iterations and time consumed
iter = temp.iter;
itrej = temp.itrej; 
ns = temp.ns;
pcgit = temp.pcgit; 
lpit = temp.lpit;
totaltime = array2table(sum(temp.time{:,:},1));
totaltime.Properties.VariableNames = temp.time.Properties.VariableNames;
temp.time = [temp.time; totaltime]; 
temp.rownames = [temp.rownames; "Total"]; 
% ---------- %

% Setting the table of times %
time = varfun(@(x) num2str(x, ['%' sprintf('.%df', 2)]), temp.time);
time.Properties.VariableNames = {'PreFilter','Filter','SetupK','SetupPrec','Systems','Grad','LP','Fix','Proj','Other','Total'};
time.Properties.RowNames = temp.rownames;
% ---------- %

end