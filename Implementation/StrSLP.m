function [xstar,xstarDsgn,u,opstop,iter,itrej,F,vol,ns,pcgit,lpit,time,xhist] = StrSLP(str,volfrac,p,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,strDsgn,x0,xfil0,W,w,k,kv,iK,jK,perm,Dofs,time)
% StrSLP solves the topology optimization problem of minimum compliance for structures using a Sequential Linear Programming (SLP) algorithm. %
% INPUT: str - structure with the problem data.
%        volfrac - maximum volume fraction of the domain that the structure can occupy.
%        p - penalty parameter for the SIMP model. 
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
%        genK - explicit generation of the global stiffness matrix (true = yes, false = no).
%        mr - structure that contains parameters used by the multiresolution method.
%        strDens - structure with density mesh data for multiresolution. 
%        strDsgn - structure with design variable mesh data for multiresolution.
%        x0 - initial design variables.
%        xfil0 - initial densities vector.
%        W - matrix containing the weight factors for each finite element, associated to the filter. 
%        w - vector with the sum of the weight factors for each finite element. 
%        k - element stiffness matrix. 
%        kv - element stiffness matrix converted into a vector. 
%        iK - row indexes with nonzero elements of the global stiffness matrix K.
%        jK - column indexes with nonzero elements of the global stiffness matrix K.
%        perm - approximate minimum degree permutation vector of K.
%        Dofs - matrix with the indexes of the degrees of freedom for each element.
%        time - structure with the time consumed on each part of the program.
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
% ---------- %

% Some initial values %
x = x0; % Initial design variables
delta = slp.delta0;
V = str.l*str.h*str.w; % Domain's volume
if(mr.op) % Multiresolution
    nelem = strDsgn.nelem; 
    velem = (V/strDens.nelem)*ones(strDens.nelem,1); % Element volumes vector
else
    nelem = str.nelem;    
    velem = (V/nelem)*ones(nelem,1); % Element volumes vector
end 
Vmax = volfrac*V; % Maximum volume that the structure can occupy
theta1 = 1; % Penalty parameter for the merit function
theta2 = 1; % Penalty parameter for the merit function
thetaMax = 1; % Penalty parameter for the merit function
countS = 0; % Number of times the step norm is less or equal than 'tolS'
countG = 0; % Number of times the norm of the projected gradient is less or equal than 'tolG'
countF = 0; % Number of times the actual reduction of the objective function is less or equal than 'tolF'
iter = 0; % Number of iterations
itrej = 0; % Number of rejected steps
pcgit = zeros(slp.maxiter,1); % Number of iterations of PCG for each system solved
lpit = zeros(slp.maxiter,1); % Number of iterations of linprog for each LP solved
xfil = xfil0; % Initial densities vector
if(genK)
    xhist = cell(slp.maxiter+1,1);
    %xhist{1} = xfil;
else
    xhist = [];
end
str.freeDensG = str.freeDens; % Elements that will be in the computation of the gradient 
% ---------- % 

% Finite element analysis %
tic;
if (genK)
    if(opsolver == 3 || opsolver == 4) % Using multigrid
        Ac = GlobalStiffMatrix(str,kv,iK,jK,xfil,p,elem,emin,mr,strDens);
    else
        K = GlobalStiffMatrix(str,kv,iK,jK,xfil,p,elem,emin,mr,strDens);
        if(elem.fix ~= 0) % Recalculate perm if the size of K changed when fixing dofs
            perm = amd(K); 
        end
    end
end
time.setupK = time.setupK + toc; 
u = zeros(3*str.nnodes,1);
if(opsolver == 3) % Using AMG 
    tic; 
    if (genK)
        [Ac,P,M,L,U,R,perm,mg.ngrids] = SetupAMG(Ac,mg);
    else
        [An,P,M,R,perm,mg.ngrids] = SetupAMGWOK(str,mg,kv,xfil,p,emin,elem,mr);
    end
    time.setupPrec = time.setupPrec + toc; 
    tic;
    if (genK)
        [u(str.freedofs),pcgit(1)] = SolveMG(u(str.freedofs),str.f,Ac,P,M,L,U,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,pcgp);
    else 
        [u(str.freedofs),pcgit(1)] = SolveMGWOK(u(str.freedofs),An,P,M,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,kv,xfil,p,emin,pcgp,str,elem,mr);
    end
    time.system = time.system + toc;
elseif(opsolver == 4) % Using GMG
    tic;
    if (genK)
        [Ac,P,M,L,U,R,perm] = SetupGMG(Ac,str,mg,elem);
    else
        [An,P,M,R,perm] = SetupGMGWOK(str,mg,kv,xfil,p,emin,elem,mr);
    end
    time.setupPrec = time.setupPrec + toc;
    tic;
    if (genK)
        [u(str.freedofs),pcgit(1)] = SolveMG(u(str.freedofs),str.f,Ac,P,M,L,U,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,pcgp);
    else
        [u(str.freedofs),pcgit(1)] = SolveMGWOK(u(str.freedofs),An,P,M,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,kv,xfil,p,emin,pcgp,str,elem,mr);
    end
    time.system = time.system + toc; 
elseif ((opsolver == 1) && (~genK))
    [u(str.freedofs),pcgit(1),tp,ts] = NodalDispWOK(str,kv,xfil,p,emin,u(str.freedofs),pcgp,elem,mr); % Initial nodal displacements vector
    time.setupPrec = time.setupPrec + tp;
    time.system = time.system + ts;
else
    [u(str.freedofs),pcgit(1),tp,ts] = NodalDisp(K,str.f,perm,opsolver,u(str.freedofs),pcgp); % Initial nodal displacements vector
    time.setupPrec = time.setupPrec + tp;
    time.system = time.system + ts;
end
ns = 1; % Number of linear systems solved 
% ---------- % 

% Objective function and its gradient vector %
if(mr.op && mr.interp) % Using the interpolated displacements on the density grid to calculate the objective function
    strT = strDens; 
    strT.supp = str.supp; 
    strT.loads = str.loads;
    strT = ConvertDofs(strT,elem); 
    DofsDens = DegreesOfFreedom(strDens,elem);
    ud = zeros(3*strT.nnodes,1);
    nely = strDens.nely; 
    nelxy = strDens.nelx*strDens.nely;
    for i = 1:str.nelem
        vb = zeros(mr.n^2, 1); 
        ind = 1; 
        for n1 = 0:(mr.n-1) 
            for n2 = 0:(mr.n-1)
                vb(ind) = strDens.fde(i) + n1*nelxy + n2*nely; 
                ind = ind + 1; 
            end
        end     
        v = zeros(mr.n^2, mr.n);
        for n3 = 0:(mr.n-1)
            v(:,n3+1) = vb+n3; 
        end
        v = v(:);         
        v = sort(v); % Indexes of the density elements inside the finite element i   
        uel = u(Dofs(i,:)); % Nodal displacements of the finite element i 
        ueli = mr.Idisp*uel; % Interpolated displacements of all density elements inside the finite element i
        for i2 = 1:length(v)       
            ud(DofsDens(v(i2),:)) = ueli(mr.Ldofs(i2,:)); % Interpolated displacements of a single density element            
        end
    end
    Fold = strT.f'*ud(strT.freedofs);
else
    Fold = str.f'*u(str.freedofs);
end
tic;
if (genK)
    Fgrad = Grad(str,u,p,emin,k,xfil,Dofs,W,w,opfilter,mr,strDens);
else
    Fgrad = GradWOK(str,u,p,emin,k,xfil,rmin,w,opfilter,elem,mr,strDens,strDsgn);
end
time.grad = time.grad + toc; 
% ---------- % 

% Choosing variables to suppress from the problem, having fixed values % 
if(elem.fix ~= 0)
    tic;
    [elem2fix0,elem2fix1,dofs2fix,dsgn2fix,dens2fix] = Elem2Fix(str,strDsgn,strDens,x,xfil,mr,Dofs,genK,elem,Fgrad);

    strDsgn.oldfixedDens = strDsgn.fixedDens; % Save the original fixed design variables 
    strDsgn.fixedDens = union(strDsgn.oldfixedDens,dsgn2fix);
    strDsgn.freeDens = setdiff(1:strDsgn.nelem,strDsgn.fixedDens);
    if(mr.op && mr.n ~= mr.d)
        strDens.oldfixedDens = strDens.fixedDens; % Save the original fixed densities 
        strDens.fixedDens = union(strDens.oldfixedDens,dens2fix);
        strDens.freeDens = setdiff(1:strDens.nelem,strDens.fixedDens);
    end
    str.oldfixedDens = str.fixedDens; % Save the original fixed displacement elements
    str.fixedDens = union(str.oldfixedDens,[elem2fix0,elem2fix1]);
    str.freeDens = setdiff(1:str.nelem,str.fixedDens);
    str.freeDensG = setdiff(1:str.nelem,union(str.oldfixedDens,elem2fix0)); % Elements that will remain in the computation of the gradient 
    
    if(elem.fixdofs ~= 0)
        str.freedofs = setdiff(1:3*str.nnodes,union(str.suppindex,dofs2fix)); 
        str.f = str.oldf(str.freedofs);
        u(dofs2fix) = 0;
        resetup = true;
    else
        resetup = false;
    end
    
    time.fix = time.fix + toc; 
else
    resetup = false;
end
% ---------- % 

% Assigning fixed design variables % 
if(mr.op)
    fixedDsgn = strDsgn.fixedDens; 
    freeDsgn = strDsgn.freeDens;
else
    fixedDsgn = str.fixedDens;
    freeDsgn = str.freeDens;
end
% ---------- %

% LP constraints % 
if(opfilter == 0)
    A = velem'/Vmax; 
elseif(opfilter == 1 || opfilter == 2)
    if (genK)
        A = (ApplyFilter(velem,W,w,opfilter,true))'/Vmax;
    else
        A = ApplyFilterWOK(velem,str,rmin,w,opfilter,true,mr,strDens,strDsgn)'/Vmax;
    end
end
b = (velem'*xfil - Vmax)/Vmax; 
% ---------- % 

% Lower and upper bounds for the LP problem variables %
sL = max(-delta, -x);
sU = min(delta, 1-x);
sL(fixedDsgn) = 0.0;
sU(fixedDsgn) = 0.0;
% ---------- % 

% Choosing the linprog parameter and some options according to the LP solver %
if (slp.lpsolver == 1) 
    LPoptions = optimset('Display','off','MaxIter',100000,'Algorithm','interior-point','TolCon',1e-8);
else
    LPoptions = optimset('Display','off','MaxIter',100000,'Algorithm','dual-simplex'); 
end
% ---------- %

% Checking if the initial solution is already optimal.
Lgrad = Fgrad(freeDsgn);
gradpbox = min(1,max(0.0,x(freeDsgn)-Lgrad)) - x(freeDsgn); 
gpnorm = norm(gradpbox,inf);
if ((gpnorm<1e-12)&&(b<1e-12)&&(max(x-1.0)<1e-12)&&(min(x)>-1e-12))
    % The KKT conditions were approximately satisfied by the initial solution %
    countG = slp.maxcount;
    countF = slp.maxcount;
    F = Fold;
    lambda.ineqlin = 0.0;
end
% ---------- %

% SLP algorithm %
while((countG < slp.maxcount || countF < slp.maxcount) && (countS < slp.maxcount) && (iter < slp.maxiter))
    it = iter + 1; % Current iteration

    % Solving the LP problem %
    tic;
    if(volineq)
        [s,fval,flag,outLP,lambda] = linprog(Fgrad,A,-b,[],[],sL,sU,LPoptions);
    else
        [s,fval,flag,outLP,lambda] = linprog(Fgrad,[],[],A,-b,sL,sU,LPoptions);
    end
    if(flag ~= 1)
        A2 = [A -1 1];
        c = [zeros(nelem,1); 1; 1];
        sL2 = [max(-0.8*delta, -x); 0; 0];
        sU2 = [min(0.8*delta, 1-x); inf; inf];
        sL2(fixedDsgn) = 0.0; 
        sU2(fixedDsgn) = 0.0;
        if(volineq)
            [s,~,~,outLP,lambda] = linprog(c,A2,-b,[],[],sL2,sU2,LPoptions);
        else
            [s,~,~,outLP,lambda] = linprog(c,[],[],A2,-b,sL2,sU2,LPoptions);
        end
        s = s(1:nelem);
        fval = Fgrad'*s; 
    elseif (outLP.iterations==0)
        s = zeros(nelem,1);
    end
    time.LP = time.LP + toc;
    lpit(ns) = outLP.iterations;

    % Step norm %
    snorm = norm(s,inf);

    % Saving some values and updating x %
    xold = x;
    x = x + s;
    bold = b;
    deltaold = delta;

    % Applying the filter %
    tic;
    if (genK)
        xfil = ApplyFilter(x,W,w,opfilter,false);
        %xhist{iter+1} = xfil;
    else
        xfil = ApplyFilterWOK(x,str,rmin,w,opfilter,false,mr,strDens,strDsgn);
    end
    time.filter = time.filter + toc;

    % Updating b %
    b = (velem'*xfil - Vmax)/Vmax;

    % Finite element analysis (solving the linear system) %  
    if (genK)
        tic;
        if(opsolver == 3 || opsolver == 4)
            Ac{1} = GlobalStiffMatrix(str,kv,iK,jK,xfil,p,elem,emin,mr,strDens);
        else
            K = GlobalStiffMatrix(str,kv,iK,jK,xfil,p,elem,emin,mr,strDens);
            % Recalculate perm if the size of K changed when fixing dofs
            if(resetup)
                perm = amd(K); 
            end
        end
        time.setupK = time.setupK + toc;
    end 
    if(pcgp.opu0)
        uold = u(str.freedofs);
    else
        uold = zeros(size(u(str.freedofs)));
    end
    if(opsolver == 3 || opsolver == 4) % Using Multigrid
        if (genK)
            tic;
            if(resetup) % Setup multigrid again if we choosed displacements to be fixed
                if(opsolver == 4)
                    [Ac,P,M,L,U,R,perm] = SetupGMG(Ac{1},str,mg,elem);
                elseif(opsolver == 3)
                    [Ac,P,M,L,U,R,perm,mg.ngrids] = SetupAMG(Ac{1},mg);
                end
                resetup = false; 
            else
                for j = 1:(mg.ngrids-1) 
                    Ac{j+1} = P{j}'*Ac{j}*P{j}; 
                    if(mg.smoother == 0) % Jacobi
                        M{j} = (1/mg.omega)*diag(diag(Ac{j}));
                    elseif(mg.smoother == 1) % Gauss-Seidel/SOR
                        M{j} = (1/mg.omega)*diag(diag(Ac{j})) + tril(Ac{j},-1);
                    elseif(mg.smoother == 2) % SSOR 
                        L{j} = (1/mg.omega)*diag(diag(Ac{j})) + tril(Ac{j},-1); 
                        U{j} = diag(diag(Ac{j}).^-1)*L{j}';
                        L{j} = (mg.omega/(2-mg.omega))*L{j}; 
                    end
                end        
                if (mg.opchol == 1)
                    if (mg.opcond == 1)
                        R = ichol(Ac{mg.ngrids}(perm,perm), struct('diagcomp',0.1));
                    else
                        R = diag(diag(Ac{mg.ngrids}));
                    end
                else
                    R = chol(Ac{mg.ngrids}(perm,perm));
                end
            end
            time.setupPrec = time.setupPrec + toc;
            tic; 
            [u(str.freedofs),pcgit(ns+1)] = SolveMG(uold,str.f,Ac,P,M,L,U,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,pcgp);
            time.system = time.system + toc;
        else
            tic;
            [An,M] = GlobalStiffMatrixWOK(P,str,kv,xfil,p,emin,mg.omega,mg.ngrids,mg.Anbatch,elem,mr);
            if (mg.opchol == 1)
                if (mg.opcond == 1)        
                    R = ichol(An(perm,perm), struct('diagcomp',0.1));
                else
                    R = diag(diag(An));
                end
            else
                R = chol(An(perm,perm));
            end            
            time.setupPrec = time.setupPrec + toc;
            tic; 
            [u(str.freedofs),pcgit(ns+1)] = SolveMGWOK(uold,An,P,M,R,perm,mg.ngrids,mg.cycle,mg.smoother,mg.smoothiter1,mg.smoothiter2,mg.tolMG,mg.maxiterMG,mg.opchol,mg.opcond,kv,xfil,p,emin,pcgp,str,elem,mr);
            time.system = time.system + toc;
        end
    elseif ((opsolver == 1) && (~genK))
        [u(str.freedofs),pcgit(ns+1),tp,ts] = NodalDispWOK(str,kv,xfil,p,emin,uold,pcgp,elem,mr);
        time.setupPrec = time.setupPrec + tp;
        time.system = time.system + ts;
    else
        [u(str.freedofs),pcgit(ns+1),tp,ts] = NodalDisp(K,str.f,perm,opsolver,uold,pcgp);
        time.setupPrec = time.setupPrec + tp;
        time.system = time.system + ts;
    end
    ns = ns+1;

    % Updating the objective function %
    if(mr.op && mr.interp) % Using the interpolated displacements on the density grid to calculate the objective function
        for i = 1:str.nelem
            vb = zeros(mr.n^2, 1); 
            ind = 1; 
            for n1 = 0:(mr.n-1) 
                for n2 = 0:(mr.n-1)
                    vb(ind) = strDens.fde(i) + n1*nelxy + n2*nely; 
                    ind = ind + 1; 
                end
            end     
            v = zeros(mr.n^2, mr.n);
            for n3 = 0:(mr.n-1)
                v(:,n3+1) = vb+n3; 
            end
            v = v(:);         
            v = sort(v); % Indexes of the density elements inside the finite element i   
            uel = u(Dofs(i,:)); % Nodal displacements of the finite element i
            ueli = mr.Idisp*uel; % Interpolated displacements of all density elements inside the finite element i
            for i2 = 1:length(v)   
                ud(DofsDens(v(i2),:)) = ueli(mr.Ldofs(i2,:)); % Interpolated displacements of a single density element          
            end
        end
        F = strT.f'*ud(strT.freedofs);
    else
        F = str.f'*u(str.freedofs);
    end

    % Predicted and actual reduction of the objective function and of the infeasibility % 
    predopt = -fval;
    predfsb = abs(bold) - abs(A*s + bold);
    aredopt = Fold - F;
    aredfsb = abs(bold) - abs(b);

    % Updating the penalty parameter (theta) for the merit function %
    thetaMin = min(theta1,theta2);
    thetaLarge = (1 + (1e6/((iter+1)^(1.1))))*thetaMin;
    if(predopt > 0.5*predfsb)
        thetaSup = 1;
    %elseif ((predopt<0)&&(aredopt<0)&&(aredfsb>0))
    %    thetaSup = max(0.5*aredfsb/(aredfsb-aredopt),0.1*predfsb/(predfsb-predopt));
    else
        thetaSup = (0.5*predfsb)/(predfsb-predopt);
    end
    theta = min(thetaLarge,thetaSup);
    theta = min(theta,thetaMax);
    theta2 = theta1;
    theta1 = theta;

    % Predicted reduction of the merit function %
    pred = theta*predopt + (1-theta)*predfsb;

    % Actual reduction of the merit function %
    ared = theta*aredopt + (1-theta)*aredfsb;

    % Testing if the new point will be accepted or rejected %
    if(ared < slp.eta*pred) % Rejected
        itrej = itrej + 1;
        delta = max(slp.alphaR*snorm, 0.1*deltaold);
        thetaMax = theta;
        x = xold;
        b = bold;
        F = Fold;
    else % Accepted
        if (ared >= slp.rho*pred) 
            delta = min(slp.alphaA*deltaold, 1.0);
        end
        delta = max(delta, 0.0001);

        thetaMax = 1;
        Fold = F;

        % Updating the objective function gradient vector %
        tic;
        if (genK)
            Fgrad = Grad(str,u,p,emin,k,xfil,Dofs,W,w,opfilter,mr,strDens);
        else
            Fgrad = GradWOK(str,u,p,emin,k,xfil,rmin,w,opfilter,elem,mr,strDens,strDsgn);
        end
        time.grad = time.grad + toc;

        % Lagrangian gradient projected onto the box %
        if(volineq)
            Lgrad = Fgrad(freeDsgn) + lambda.ineqlin*A(freeDsgn)';
        else
            Lgrad = Fgrad(freeDsgn) + lambda.eqlin*A(freeDsgn)';
        end
        gradpbox = min(1,max(0.0,x(freeDsgn)-Lgrad)) - x(freeDsgn);
        gpnorm = norm(gradpbox,inf);

        iter = iter + 1;
        
        % Update the fixed variables %   
        if(elem.fix == 2 || (elem.fix == 3 && mod(iter,elem.fixit) == 0) || (elem.fix == 4 && iter == elem.fixit))
            tic;            
            [elem2fix0,elem2fix1,dofs2fix,dsgn2fix,dens2fix] = Elem2Fix(str,strDsgn,strDens,x,xfil,mr,Dofs,genK,elem,Fgrad);

            % Fixing some variables
            strDsgn.fixedDens = union(strDsgn.oldfixedDens,dsgn2fix);
            strDsgn.freeDens = setdiff(1:strDsgn.nelem,strDsgn.fixedDens);
            if(mr.n ~= mr.d)
                strDens.fixedDens = union(strDens.oldfixedDens,dens2fix);
                strDens.freeDens = setdiff(1:strDens.nelem,strDens.fixedDens);
            end
            str.fixedDens = union(str.oldfixedDens,[elem2fix0,elem2fix1]);
            str.freeDens = setdiff(1:str.nelem,str.fixedDens);
            str.freeDensG = setdiff(1:str.nelem,union(str.oldfixedDens,elem2fix0)); % Elements that will remain in the computation of the gradient 

            if(mr.op)
                fixedDsgn = strDsgn.fixedDens; 
                freeDsgn = strDsgn.freeDens;
            else
                fixedDsgn = str.fixedDens;
                freeDsgn = str.freeDens;
            end

            % Fixing degrees of freedom (to eliminate from the linear system)
            if(elem.fixdofs ~= 0)
                str.freedofs = setdiff(1:3*str.nnodes,union(str.suppindex,dofs2fix)); 
                str.f = str.oldf(str.freedofs); 
                u(dofs2fix) = 0; 
                resetup = true;
            end 
            time.fix = time.fix + toc;
        end

        % Counting the number of times the projected gradient norm and the actual reduction of the objective function reached the tolerance value %
        if (gpnorm <= slp.tolG)
            countG = countG + 1;
        else
            countG = 0;
        end
        if (ared <= slp.tolF)
            countF = countF + 1;
        else
            countF = 0;
        end
    end
    
    % Counting the number of times the step norm reached the tolerance value %
    if (snorm <= slp.tolS) 
        countS = countS + 1;
    else
        countS = 0;
    end

    % Updating the lower and upper bounds for the LP problem variables %
    sL = max(-delta, -x);
    sU = min(delta, 1-x);
    sL(fixedDsgn) = 0.0;
    sU(fixedDsgn) = 0.0;
    
    % Display some results at each iteration % 
    fprintf('It: %3d | F: %+e | Ared: %+e | Pred: %+e | snorm: %f | delta: %f | gpnorm: %f | flagLP: %d | itLP: %3d | pcgit: %4d \n', it, F, ared, pred, snorm, deltaold, gpnorm, flag, lpit(ns-1), pcgit(ns));
end
% ---------- %

if (countG >= slp.maxcount && countF >= slp.maxcount)
    opstop = 0; % The KKT conditions were approximately satisfied and the actual reduction of the objective function was sufficiently small 
elseif (countS >= slp.maxcount)
    opstop = 1; % The step norm is sufficiently small
else
    opstop = 2; % The maximum number of iterations was reached
end

pcgit = pcgit(1:ns);
lpit = lpit(1:ns-1);

% Applying the projection to round the densities to 0 or 1
if (prj.op) 
    tic;
    if (prj.heavi)
        vcur = velem'*xfil;
        [xfil, prj] = HeavisideProj(xfil,prj,velem,vcur,vcur,false);
    end

    if ((mr.op&&(mr.n~=mr.d))||(prj.heavi))
        if (mr.n~=mr.d)
            opf = 0;
        else
            opf = opfilter;
        end
        if (genK)
            Fgrad = Grad(str,u,p,emin,k,xfil,Dofs,W,w,opf,mr,strDens);
        else
            Fgrad = GradWOK(str,u,p,emin,k,xfil,rmin,w,opfilter,elem,mr,strDens,strDsgn);
        end
        A = velem'/Vmax; 
    end    
    if(volineq)
        Lgrad = Fgrad + lambda.ineqlin*A';
    else
        Lgrad = Fgrad + lambda.eqlin*A';
    end
    xstar = GradProj(xfil,Lgrad,prj,Vmax,velem);
    time.proj = time.proj + toc;
else
    xstar = xfil; 
end

xstarDsgn = x; 
vol = (velem'*xstar)/V;

% Displaying information about the convergence criteria % 
fprintf('StrSLP summary\n')
fprintf('  opstop: %1d | countG: %2d | countF: %2d | countS: %2d | iter: %3d | it rej: %3d | F: %e | gpnorm: %f | lin sys: %5d | pcg it: %5d | lp it: %5d \n', opstop, countG, countF, countS, iter, itrej, F, gpnorm, sum(ns), sum(pcgit), sum(lpit));

end