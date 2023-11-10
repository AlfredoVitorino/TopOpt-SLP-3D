function [strDens,strDsgn,x0] = StartMR(str,volfrac,p,p123,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr)
% StartMR creates the density and design variable meshes data and gets the initial design variables vector for multiresolution. %
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
% OUTPUT: strDens - structure with density elements grid data for multiresolution. 
%         strDsgn - structure with design variable elements grid data for multiresolution.
%         x0 - initial design variables.
% ---------- %

strDens.l = str.l; 
strDens.h = str.h; 
strDens.w = str.w; 
strDens.nelx = str.nelx*(mr.n); 
strDens.nely = str.nely*(mr.n); 
strDens.nelz = str.nelz*(mr.n);
strDens.el = strDens.l/strDens.nelx; 
strDens.eh = strDens.h/strDens.nely; 
strDens.ew = strDens.w/strDens.nelz; 
strDens.nx = strDens.nelx+1; 
strDens.ny = strDens.nely+1;
strDens.nz = strDens.nelz+1;
strDens.nelem = strDens.nelx*strDens.nely*strDens.nelz;
strDens.nnodes = strDens.nx*strDens.ny*strDens.nz; 
strDens.E = str.E; 
strDens.nu = str.nu;
strDsgn.l = str.l; 
strDsgn.h = str.h; 
strDsgn.w = str.w; 
strDsgn.nelx = str.nelx*(mr.d); 
strDsgn.nely = str.nely*(mr.d); 
strDsgn.nelz = str.nelz*(mr.d);
strDsgn.el = strDsgn.l/strDsgn.nelx; 
strDsgn.eh = strDsgn.h/strDsgn.nely; 
strDsgn.ew = strDsgn.w/strDsgn.nelz; 
strDsgn.nx = strDsgn.nelx+1; 
strDsgn.ny = strDsgn.nely+1;
strDsgn.nz = strDsgn.nelz+1;
strDsgn.nelem = strDsgn.nelx*strDsgn.nely*strDsgn.nelz;
strDsgn.nnodes = strDsgn.nx*strDsgn.ny*strDsgn.nz; 
strDsgn.E = str.E; 
strDsgn.nu = str.nu;

% Calculating the indexes of the first design variable (fdv) and the first density element (fde) inside each finite element %
strDsgn.fdv = zeros(str.nelem,1);
strDens.fde = zeros(str.nelem,1);
nelyn = strDens.nely; 
nelxyn = strDens.nelx*strDens.nely;
nelyd = strDsgn.nely; 
nelxyd = strDsgn.nely*strDsgn.nelx;
jn = 1;
jd = 1;
for i = 1:str.nelem
    strDens.fde(i) = jn;
    jn = jn + mr.n;
    if (mod(jn-1,nelyn) == 0)
        jn = jn+(mr.n-1)*nelyn;
        if (mod(jn,nelxyn) == 1)
            jn = jn+(mr.n-1)*nelxyn;
        end
    end
    strDsgn.fdv(i) = jd;
    jd = jd + mr.d;
    if (mod(jd-1,nelyd) == 0)
        jd = jd+(mr.d-1)*nelyd;
        if (mod(jd,nelxyd) == 1)
            jd = jd+(mr.d-1)*nelxyd;
        end
    end
end
% ---------- %

% Elements with fixed design variables or densities (void or solid regions) %
if(~isempty(str.fixedDens))
    fixDens = zeros(strDsgn.nelem,1);
    fixVal = zeros(strDsgn.nelem,1);
    for i = 1:size(str.fix,1)   
        ex0 = (str.fix(i,1)-1)*(mr.d) + 1; % Initial element layer of the fixed design variables in the x direction 
        exf = str.fix(i,2)*(mr.d); % Final element layer of the fixed design variables in the x direction
        ey0 = (str.fix(i,3)-1)*(mr.d) + 1; % Initial element layer of the fixed design variables in the y direction
        eyf = str.fix(i,4)*(mr.d); % Final element layer of the fixed design variables in the y direction
        ez0 = (str.fix(i,5)-1)*(mr.d) + 1; % Initial element layer of the fixed design variables in the z direction
        ezf = str.fix(i,6)*(mr.d); % Final element layer of the fixed design variables in the z direction
        for elz = ez0:ezf
            for elx = ex0:exf
                for ely = ey0:eyf
                    el = ely + (elx-1)*strDsgn.nely + (elz-1)*strDsgn.nelx*strDsgn.nely;
                    fixDens(el) = 1;
                    if(str.fix(i,7) == 1) % solid region
                        fixVal(el) = 1.0;
                    end
                end
            end
        end
    end
    dens = 1:strDsgn.nelem;
    strDsgn.freeDens = dens(fixDens == 0); % elements with design variables not fixed
    strDsgn.fixedDens = dens(fixDens == 1); % elements with design variables fixed 
    strDsgn.fixedDensVal = fixVal(strDsgn.fixedDens); % design variable values of the elements with fixed design variables
    
    if(mr.n ~= mr.d) 
        fixDens = zeros(strDens.nelem,1);
        fixVal = zeros(strDens.nelem,1);
        for i = 1:size(str.fix,1)   
            ex0 = (str.fix(i,1)-1)*(mr.n) + 1; % Initial element layer of the fixed densities in the x direction 
            exf = str.fix(i,2)*(mr.n); % Final element layer of the fixed densities in the x direction
            ey0 = (str.fix(i,3)-1)*(mr.n) + 1; % Initial element layer of the fixed densities in the y direction
            eyf = str.fix(i,4)*(mr.n); % Final element layer of the fixed densities in the y direction
            ez0 = (str.fix(i,5)-1)*(mr.n) + 1; % Initial element layer of the fixed densities in the z direction
            ezf = str.fix(i,6)*(mr.n); % Final element layer of the fixed densities in the z direction
            for elz = ez0:ezf
                for elx = ex0:exf
                    for ely = ey0:eyf
                        el = ely + (elx-1)*strDens.nely + (elz-1)*strDens.nelx*strDens.nely;
                        fixDens(el) = 1;
                        if(str.fix(i,7) == 1) % solid region
                            fixVal(el) = 1.0;
                        end
                    end
                end
            end
        end
        dens = 1:strDens.nelem;
        strDens.freeDens = dens(fixDens == 0); % elements with densities not fixed
        strDens.fixedDens = dens(fixDens == 1); % elements with densities fixed 
        strDens.fixedDensVal = fixVal(strDens.fixedDens); % density values of the elements with fixed densities
    else
        strDens.freeDens = strDsgn.freeDens;
        strDens.fixedDens = strDsgn.fixedDens; 
        strDens.fixedDensVal = strDsgn.fixedDensVal;
    end
else
    strDsgn.freeDens = 1:strDsgn.nelem;
    strDsgn.fixedDens = []; 
    strDsgn.fixedDensVal = [];
    strDens.freeDens = 1:strDens.nelem;
    strDens.fixedDens = []; 
    strDens.fixedDensVal = []; 
end
% ---------- %

if(mr.x0) % Solve the problem on the coarse grid and use the solution as initial design variables for the problem on the fine grid  
    mr.op = false;
    ts = slp.tolS;
    tg = slp.tolG;
    slp.tolS = 100*slp.tolS;
    slp.tolG = 100*slp.tolG;
    disp('Solving the problem on the coarse grid to obtain the initial design variables for the problem on the fine grid.');
    [~,xstarDsgn,~,~,~,~,~,~,~,~,~,~,~,~] = StartStrSLP(str,volfrac,p,p123,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr);    
    slp.tolS = ts;                          
    slp.tolG = tg;  
    x0 = zeros(strDsgn.nelem,1);    
    nely = strDsgn.nely; 
    nelxy = strDsgn.nely*strDsgn.nelx;
    j = 1;
    for i = 1:str.nelem 
        v = zeros(mr.d^3, 1);
        ind = 1;     
        for n1 = 0:(mr.d-1) 
            for n2 = 0:(mr.d-1)
                for n3 = 0:(mr.d-1)
                    indv = j + n1*nelxy + n2*nely + n3;
                    v(ind) = indv;
                    ind = ind+1;
                end
            end
        end 
        x0(v) = xstarDsgn(i); 
        j = j+mr.d;
        if (mod(j-1,nely) == 0)
            j = j+(mr.d-1)*nely;
            if (mod(j,nelxy) == 1)
                j = j+(mr.d-1)*nelxy;
            end
        end
    end 
    mr.op = true; 
    disp('Solving the problem on the fine grid.');
else
    x0 = volfrac*ones(strDsgn.nelem,1);
    if(~isempty(strDsgn.fixedDens)) % Assign fixed design variables
        x0(strDsgn.fixedDens) = strDsgn.fixedDensVal;
        velemDsgn = ((str.l*str.h*str.w)/strDsgn.nelem)*ones(strDsgn.nelem,1); % Volume of the elements on the design variables grid
        Vfix = velemDsgn(strDsgn.fixedDens)'*x0(strDsgn.fixedDens); % Volume occupied by fixed elements
        x0(strDsgn.freeDens) = (volfrac*(str.l*str.h*str.w)-Vfix)/(sum(velemDsgn(strDsgn.freeDens))); % Redistributing initial design variables of the free elements
    end
end

end