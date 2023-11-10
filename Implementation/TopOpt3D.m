% Main program to solve the three-dimensional structural topology optimization problem of minimum compliance, using a sequential linear programming algorithm %

% Problem data (str)
% load('Tests/cb30x10x10half.mat'); % cantilever beam
% load('Tests/cb48x16x16.mat');
% load('Tests/cb48x16x16half.mat');
% load('Tests/cb60x20x20.mat'); 
% load('Tests/cb72x24x24.mat'); 
% load('Tests/cb90x30x30.mat'); 
% load('Tests/cb96x32x32.mat'); 
% load('Tests/cb96x32x32half.mat');
% load('Tests/cb120x40x40.mat');
% load('Tests/cb120x40x40half.mat');
% load('Tests/cb144x48x48.mat');
% load('Tests/cb144x48x48half.mat');
% load('Tests/cb150x50x50.mat');
% load('Tests/cb168x56x56.mat');
% load('Tests/cb180x60x60.mat');
% load('Tests/cb192x64x64.mat');
% load('Tests/cb192x64x64half.mat');
% load('Tests/cb60x30x2.mat');
load('Tests/mbb96x16x16quarter.mat'); % mbb beam
% load('Tests/mbb120x20x20quarter.mat');
% load('Tests/mbb192x32x32quarter.mat');
% load('Tests/mbb240x40x40quarter.mat');
% load('Tests/mbb288x48x48quarter.mat');
% load('Tests/mbb384x64x64quarter.mat');
% load('Tests/mbb96x16x1half.mat');
% load('Tests/tt16x80x16quarter.mat'); % transmission tower
% load('Tests/tt32x160x32quarter.mat');
% load('Tests/tt48x240x48quarter.mat');
% load('Tests/tt64x320x64quarter.mat');
% load('Tests/eb60x20x30.mat'); % engine bracket
% load('Tests/eb48x16x24half.mat'); 
% load('Tests/hc24x24x24quarter.mat'); % cube with concentrated load at the bottom 
% load('Tests/hc32x16x32quarter.mat');
% load('Tests/hc96x48x96quarter.mat');
% load('Tests/bt16x64x16.mat'); % building with torsion
% load('Tests/bt32x128x32.mat');  
% load('Tests/ls48x48x16half.mat'); % L-shaped beam
% load('Tests/ls144x144x48half.mat');
% load('Tests/ls192x192x64half.mat');
% load('Tests/bd192x48x24quarter.mat'); str.f = str.f/norm(str.f); % bridge

% Maximum volume fraction of the domain that the structure can occupy
volfrac = 0.2; % test
% volfrac = 0.2; % cantilever beam 
% volfrac = 0.2; % mbb beam
% volfrac = 0.15; % transmission tower 
% volfrac = 0.15; % engine bracket 
% volfrac = 0.16; % 1/2 cube with concentrated load at the bottom 
% volfrac = 0.1; % building with torsion
% volfrac = 0.15; % bridge (includes the fixed elements)
% volfrac = 0.1; % L-shaped beam (includes the fixed elements)

% Filter radius
rmin = 1.5; % test
% rmin = 1.4; % cb30x10x10 % cantilever beam
% rmin = 1.5; % cb48x16x16
% rmin = 1.5; % cb60x20x20
% rmin = 1.5; % cb72x24x24
% rmin = 2.0; % cb90x30x30
% rmin = 2.2; % cb96x32x32 
% rmin = 2.5; % cb120x40x40 
% rmin = 3.5; % cb144x48x48 
% rmin = 4.0; % cb150x50x50
% rmin = 4.2; % cb168x56x56
% rmin = 4.7; % cb180x60x60 
% rmin = 5.2; % cb192x64x64
% rmin = 2.0; % mbb96x16x16 % mbb beam
% rmin = 2.5; % mbb120x20x20
% rmin = 5.0; % mbb192x32x32
% rmin = 5.5; % mbb240x40x40 
% rmin = 6.0; % mbb288x48x48
% rmin = 8.0; % mbb384x64x64
% rmin = 1.1; % tt16x80x16 % transmission tower
% rmin = 2.0; % tt32x160x32
% rmin = 2.2; % tt40x200x40
% rmin = 1.5; % eb60x20x30 % engine bracket
% rmin = 1.5; % eb48x16x24
% rmin = 1.5; % hc32x16x32 % half cube with concentrated load at the bottom
% rmin = 3.0; % hc96x48x96
% rmin = 1.5; % bt16x64x16 % building with torsion
% rmin = 3.0; % bt32x128x32 
% rmin = 1.5; % ls48x48x16 % L-shaped beam
% rmin = 3.0; % ls144x144x48
% rmin = 1.5; % bd192x48x24 % bridge

p = 3; % Penalty parameter for the SIMP model
p123 = false; % Defines if we will set p = 1, 2 and 3 (p123 = true) or use only one value of p (p123 = false)
femin = 1.0e-6; % Stiffness factor of the void material (E_(void) = femin*E)
emin = femin*str.E; % Young's modulus of the void material
opfilter = 1; % Filter option (0 = no filter, 1 = weighted average density filter, 2 = average density filter) 
volineq = true; % Defines the type of the volume constraint (true = "less than or equal to", false = "equal to")

% Finite element characteristics and parameters 
elem.deg = 1; % Element degree (polynomial degree of the shape functions). It is possible to use degree 1, 2 or 3
elem.type = 1; % Element type (1 = Lagrange, 2 = Serendipity)
elem.increasedeg = false; % Defines if we will solve the problem increasing the degree of the elements until reach 'elem.maxdeg'
elem.maxdeg = 2; % Maximum element degree 
elem.fix = 4; % Defines the way we choose variables to be fixed when increasing the element degree
%(0 = don't choose variables to be fixed, 1 = choose only at the first iteration, 2 = choose in all iterations, 3 = choose every 'elem.fixit' iterations, 4 = choose in iterations 1 and 'elem.fixit')
elem.fixit = 5; % Number of outer iterations to be accomplished until choose again the variables that will be fixed
elem.fixtl = 1e-6; % Tolerance to find approximately void elements (the ones full of densities <= elem.fixtl)
elem.fixtu = 0.9; % Tolerance to find approximately solid elements (the ones full of densities >= elem.fixtu)
elem.fixnb = false; % Defines if we will choose all void elements to be fixed (fixnb = true) or just the ones surrounded by void elements (fixnb = false)
elem.fixdofs = 1; % Defines if we will fix displacements of void elements eliminating the degrees of freedom from the linear system 
%(0 = don't fix any displacements, 1 = don't fix the displacements of nodes shared by neighbor elements, 2 = fix the displacements of all nodes of void elements
elem.fixdsgn = false; % Defines if we will choose void elements looking for the ones full of approximately zero design variables (fixdsgn = true) or densities (fixdsgn = false)
elem.ltdsgn = 0.3; % Maximum value a neighbor may have when fixing design variables to zero
elem.utdsgn = 0.7; % Minimum value a neighbor may have when fixing design variables to one
elem.tolgdsgn = 1e-6; % Tolerance for the gradient when fixing design variables

% SLP parameters
slp.tolS = 1e-4; % Threshold for the step norm obtained by the SLP algorithm
slp.tolG = 1e-3; % Threshold for the projected gradient norm obtained by the SLP algorithm
slp.tolF = 5e-2; % Threshold for the actual reduction of the objective function obtained by the SLP algorithm
slp.maxcount = 3; % Number of times that the step norm (or projected gradient norm) needs to reach a value less or equal than the threshold before stopping the SLP algorithm (stopping criterion)
slp.maxiter = 500; % Maximum number of iterations (stopping criterion)
slp.maxit12 = 15;  % Maximum number of iterations for p = 1, 2
slp.lpsolver = 0; % Linear programming solver (0 = dual simplex, 1 = interior point method)
slp.delta0 = 0.1; % Initial trust region radius
slp.eta = 0.1; % Parameter used to accept the step (s is accepted if ared >= 0.1*pred)
slp.rho = 0.5; % Parameter used to increase the trust region radius (delta is increased if ared >= 0.5*pred)
slp.alphaR = 0.25; % Factor used to decrease the trust region radius when the step is rejected
slp.alphaA = 2.0; % Factor used to increase the trust region radius when the step is accepted

% Parameters used for the density projection strategy 
prj.op = false; % Defines if the densities are to be rounded to 0 or 1 at the end of the SLP algorithm.
prj.nearlim = true; % Try to round the variables to the nearest limit (0 or 1) prior to using the gradient.
prj.maxang = 89.9; % Largest angle (in degrees) allowed between (xnew-x) and -g when rounding the densities.
prj.maxit = 10; % Maximum number of calls to the SLP algorithm when prj.op = true.
prj.tolV = 0.005; % Maximum violation of the volume fraction allowed after projection.
prj.tolN = 0.01; % Maximum fraction of element density changes between consecutive projected solutions. 
prj.opfilter = 1; % Filter used after the first projection (0 = no filter, 1 = weighted average density filter, 2 = average density filter).
prj.rmin = 1.1; % Filter radius to be used when calling StrSLP after a projection. 
prj.gthresh = 1e-4; % Gradient threshold (g(i) is considered to be near 0 if abs(g(i)) <= gthresh*norm(g,inf)).
prj.ut = 0.95; % Upper limit threshold (densities greater or equal to ut are always rounded to 1)
prj.lt = 0.05; % Lower limit threshold (densities lower or equal to lt are always rounded to 0)
prj.vumin = 0.7; % Smallest density that can be rounded to 1
prj.vlmax = 0.3; % Largest density that can be rounded to 0
prj.heavi = true; % Indicates if the smooth Heaviside projection will be applied to round the densities.
prj.eta = 0.25; % Initial value of the eta parameter, used to compute the Heaviside filter function.
prj.beta = 1.0; % Initial value of the beta parameter, used to compute the Heaviside filter function.
prj.dbeta = 2.0; % Increasing factor for the beta parameter (beta(k+1) = dbeta*beta(k)).
prj.betamax = 1000.0; % Maximum value of the Heaviside beta parameter.
prj.deg = false; % Defines if we will use linear elements on the calls to the SLP algorithm after rounding the densities. 

% Parameters for considering the problem with symmetries 
sym.xy = true; % Defines if the problem domain has a symmetry with respect to the xy plane (true for cb, mbb, tt and hc problems)
sym.yz = true; % Defines if the problem domain has a symmetry with respect to the yz plane (true for mbb, tt and hc problems)
sym.xz = false; % Defines if the problem domain has a symmetry with respect to the xz plane

opsolver = 4; % Linear system solver option
%0 = Cholesky factorization, 
%1 = CG w/ preconditioner being the diagonal of K, 
%2 = CG w/ preconditioner being the incomplete Cholesky factorization of K,
%3 = CG w/ Algebraic Multigrid,
%4 = CG w/ Geometric Multigrid. 

% PCG parameters
pcgp.opu0 = true; % Adoption of the previous vector u as an initial guess for pcg
pcgp.tolPCG = 1e-8; % Convergence tolerance for pcg
pcgp.maxiterPCG = 10000; % Maximum number of iterations of pcg

genK = true; % Defines if the global stiffness matrix will be explicitly generated
% You can set genK to false whenever opsolver = 1, 3 or 4 and the problem is huge
addpath('WOK');

% Multigrid parameters 
mg.ngrids = 3; % Number of grids
mg.cycle = 1; % Cycle type (0 - Vcycle, 1 - Wcycle, 2 - FullVCycle)
mg.smoother = 0; % Smoother (0 - Jacobi, 1 - Gauss-Seidel/SOR, 2 - SSOR)
mg.omega = 0.5; % Smoother relaxation parameter 
mg.smoothiter1 = 1; % Number of pre-smooth iterations 
mg.smoothiter2 = 1; % Number of post-smooth iterations 
mg.tolMG = 1e-8; % Tolerance for CG w/ multigrid
mg.maxiterMG = 5000; % Maximum number of iterations for CG w/ multigrid 
mg.theta = 20; % Strong connections parameter 
mg.optheta = 1; % Choose if 'theta' will be used as (0 - a tolerance, 1 - the maximum number of strong connections per node)
mg.nt = 30; % Number of test space vectors
mg.tolr = 1e-1; % Tolerance for the test space generation 
mg.tolq = 1e-2; % Tolerance for the test space generation 
mg.itp = 2; % Maximum number of iterations on prolongation DPLS 
mg.kappa = 100; % Tolerance for the condition number on prolongation DPLS
mg.opcond = 1; % Preconditioner option for the test space generation and for PCG to solve the system on the coarsest grid (0 - diag, 1 - ichol)
mg.nmaxchol = 60000; % Maximum system dimension to be able to use the Cholesky factorization on the coarsest grid
mg.opchol = 0; % Action adopted when the dimension of the system on the coarsest grid is greater than 'nmaxchol' (0 - increase ngrids, 1 - use PCG)
mg.Anbatch = 4096; % Batch used in the computation of An when genK = false. An is assembled in batches of Anbatch elements.

% Multiresolution parameters
mr.op = false; % Defines if the multiresolution method will be used.
mr.n = 2; % Number of density elements on each direction of each finite element.
mr.d = 2; % Number of design variable elements on each direction of each finite element.
mr.x0 = false; % Defines if the solution of the coarse problem will be used as initial element densities.
mr.interp = false; % Defines if the displacements will be interpolated to calculate the gradient.

% Apply an SLP algorithm to solve the problem
[xstar,xstarDsgn,u,opstop,iter,itrej,F,vol,ns,pcgit,lpit,time,xhist,strDens] = StartStrSLP(str,volfrac,p,p123,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr); 

% Calculate extra results based on the solution obtained
exop.correctF = false; % Calculate the "correct" function value when multiresolution or an element degree greater than 1 is used
exop.postopt = false; % Solve the problem again on the fine mesh, with the multiresolution solution as initial guess
exop.roundsol = false; % Round the greatest densities to 1 until filling the volume and the others to 0, then calculate the objective function of the rounded solution
exop.solidsol = false; % Calculate the objective function of the full solid structure (with all densities equal to 1)
if(exop.correctF || exop.postopt || exop.roundsol || exop.solidsol)
    exres = ExtraResults(str,volfrac,p,emin,opfilter,rmin,volineq,elem,slp,opsolver,pcgp,mg,prj,genK,mr,strDens,xstarDsgn,xstar,exop);
end

% Display a 3D view of the optimal structure %
if(mr.op)
    if(sym.xy || sym.yz || sym.xz)
        xstartotal = TotalStr(strDens,xstar,sym,prj.op);
    else
        Display3D(strDens,xstar,prj.op); 
    end
else
    if(sym.xy || sym.yz || sym.xz)
        xstartotal = TotalStr(str,xstar,sym,prj.op);
    else
        Display3D(str,xstar,prj.op); 
    end
end