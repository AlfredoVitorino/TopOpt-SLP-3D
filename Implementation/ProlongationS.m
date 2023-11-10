function [P,freedofsc] = ProlongationS(nelxf,nelyf,nelzf,nelxc,nelyc,nelzc,freedofsf,elem)
% ProlongationS obtains the prolongation matrix for geometric multigrid using serendipity elements. %
% INPUT: nelxf - number of elements in the length(x) direction on the fine grid. 
%        nelyf - number of elements in the height(y) direction on the fine grid. 
%        nelzf - number of elements in the width(z) direction on the fine grid. 
%        nelxc - number of elements in the length(x) direction on the coarse grid. 
%        nelyc - number of elements in the height(y) direction on the coarse grid. 
%        nelzc - number of elements in the width(z) direction on the coarse grid.  
%        freedofsf - indexes of the free (without support) degrees of freedom for the fine nodes. 
%        elem - structure with element characteristics (deg - polynomial degree, type - element type (1 = Lagrange, 2 = serendipity).
% OUTPUT: P - prolongation matrix. 
%         freedofsc - indexes of the free (without support) degrees of freedom for the coarse nodes. 
% ---------- %

elemDeg = elem.deg; % element degree
     
ny1c = nelyc*elemDeg + 1; % number of nodes in the y-direction on layers containing vertices of the elements (on the coarse grid)
nx2c = nelxc + 1; % number of nodes in the x direction on layers not containing vertices of the elements (on the coarse grid)
ny2c = nelyc + 1; % number of nodes in the y direction on layers not containing vertices of the elements (on the coarse grid)
n3 = elemDeg - 1; % number of nodes in the interior of each edge (without counting vertices) of the element (on the coarse grid)
nxy1c = ny1c*(nelxc+1) + ny2c*n3*(nelxc); % number of nodes in an xy layer containing vertices of the element (on the coarse grid)
nxy2c = nx2c*ny2c; % number of nodes in an xy layer not containing vertices of the element (on the coarse grid)
nc = nxy1c*(nelzc+1) + nxy2c*n3*nelzc; % total number of nodes on the coarse grid

m = 8 + 12*(elemDeg-1); % total number of nodes in each element
nxf = nelxf*elemDeg+1; nyf = nelyf*elemDeg+1; nzf = nelzf*elemDeg+1; % number of node layers in each direction of the fine grid 
nf = nxf*nyf*nzf; % number of fine nodes (exceeding) 
maxind = 3*m*nf; 

% Prealocation to construct the sparse prolongation matrix P and the freedofs vector for the coarse grid
iP = zeros(maxind,1); 
jP = zeros(maxind,1);
sP = zeros(maxind,1);
freedofsc = zeros(1,3*nc);
aux = zeros(1,3*nc);

% Local coordinates of the nodes inside each element 
if(elemDeg == 1)
    localCord = [-1,0];
elseif(elemDeg == 2)
    localCord = [-1,-0.5,0,0.5];
elseif(elemDeg == 3)
    localCord = [-1,-2/3,-1/3,0,1/3,2/3];
end
localCordx = [repmat(localCord,1,nelxc),1];
localCordy = [repmat(localCord,1,nelyc),1];
localCordz = [repmat(localCord,1,nelzc),1];

% Element layers in each direction of the coarse grid 
elemLayerx = repmat(1:nelxc,2*elemDeg,1);
elemLayerx = [elemLayerx(:); nelxc];
elemLayery = repmat(1:nelyc,2*elemDeg,1);
elemLayery = [elemLayery(:); nelyc];
elemLayerz = repmat(1:nelzc,2*elemDeg,1);
elemLayerz = [elemLayerz(:); nelzc];

% Degrees of freedom of the coarse elements 
strC.nelx = nelxc; 
strC.nely = nelyc; 
strC.nelz = nelzc; 
strC.nelem = nelxc*nelyc*nelzc;
DofsC = DegreesOfFreedom(strC,elem);

% Loop in the fine nodes by layers in each direction % 
nodef = 0;
nodec = 0;
ind = 1;
for zf = 1:nzf
    for xf = 1:nxf 
        for yf = 1:nyf
            if(elemDeg == 1 || sum(mod([xf,yf,zf],elemDeg)==1) >= 2) % Only the nodes that are in the serendipity element 
                nodef = nodef + 1; % Index of the node on the fine grid (row index of P)
             
                % Find the index of the element on the coarse grid where the node is 
                elx = elemLayerx(xf);
                ely = elemLayery(yf);
                elz = elemLayerz(zf);
                elemC = (elz-1)*nelxc*nelyc + (elx-1)*nelyc + ely;  
                
                % Find the local coordinates of the node inside the element on the coarse grid 
                x1 = localCordx(xf);
                x2 = localCordy(yf);
                x3 = localCordz(zf);
                
                % Calculate the shape functions at the point (x1,x2,x3) (values of P)
                if(elemDeg == 1) % Linear element (8 shape functions)
                    phi1 = (1/8)*(1-x1)*(1-x2)*(1-x3);
                    phi2 = (1/8)*(1-x1)*(1+x2)*(1-x3);
                    phi3 = (1/8)*(1+x1)*(1-x2)*(1-x3);
                    phi4 = (1/8)*(1+x1)*(1+x2)*(1-x3);
                    phi5 = (1/8)*(1-x1)*(1-x2)*(1+x3);
                    phi6 = (1/8)*(1-x1)*(1+x2)*(1+x3);
                    phi7 = (1/8)*(1+x1)*(1-x2)*(1+x3);
                    phi8 = (1/8)*(1+x1)*(1+x2)*(1+x3);
                    v = [phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8]; 
                elseif(elemDeg == 2) % Serendipity quadratic element (20 shape functions)
                    phi1 = (1/8)*(1 - x1)*(1 - x2)*(1 - x3)*(-2 - x1 - x2 - x3); 
                    phi2 = (1/4)*(1 - x1)*(1 - x2^2)*(1 - x3);
                    phi3 = (1/8)*(1 - x1)*(1 + x2)*(1 - x3)*(-2 - x1 + x2 - x3);
                    phi4 = (1/4)*(1 - x1^2)*(1 - x2)*(1 - x3);
                    phi5 = (1/4)*(1 - x1^2)*(1 + x2)*(1 - x3);
                    phi6 = (1/8)*(1 + x1)*(1 - x2)*(1 - x3)*(-2 + x1 - x2 - x3);
                    phi7 = (1/4)*(1 + x1)*(1 - x2^2)*(1 - x3);
                    phi8 = (1/8)*(1 + x1)*(1 + x2)*(1 - x3)*(-2 + x1 + x2 - x3);
                    phi9 = (1/4)*(1 - x1)*(1 - x2)*(1 - x3^2);
                    phi10 = (1/4)*(1 - x1)*(1 + x2)*(1 - x3^2);
                    phi11 = (1/4)*(1 + x1)*(1 - x2)*(1 - x3^2);
                    phi12 = (1/4)*(1 + x1)*(1 + x2)*(1 - x3^2);
                    phi13 = (1/8)*(1 - x1)*(1 - x2)*(1 + x3)*(-2 - x1 - x2 + x3);
                    phi14 = (1/4)*(1 - x1)*(1 - x2^2)*(1 + x3);
                    phi15 = (1/8)*(1 - x1)*(1 + x2)*(1 + x3)*(-2 - x1 + x2 + x3);
                    phi16 = (1/4)*(1 - x1^2)*(1 - x2)*(1 + x3);
                    phi17 = (1/4) *(1 - x1^2)*(1 + x2)*(1 + x3);
                    phi18 = (1/8)*(1 + x1)*(1 - x2)*(1 + x3)*(-2 + x1 - x2 + x3);
                    phi19 = (1/4)*(1 + x1)*(1 - x2^2)*(1 + x3);
                    phi20 = (1/8)*(1 + x1)*(1 + x2)*(1 + x3)*(-2 + x1 + x2 + x3);
                    v = [phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9, phi10,...
                        phi11, phi12, phi13, phi14, phi15, phi16, phi17, phi18, phi19, phi20];
                elseif(elemDeg == 3) % Serendipity cubic element (32 shape functions)
                    phi1 = (1/64)*(1 - x1)*(1 - x2)*(1 - x3)*(-19 + 9*x1^2 + 9*x2^2 + 9*x3^2);
                    phi2 = (9/64)*(1 - x1)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 - x3);
                    phi3 = (9/64)*(1 - x1)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 - x3);
                    phi4 = (1/64)*(1 - x1)*(1 + x2)*(1 - x3)*(-19 + 9*x1^2 + 9*x2^2 + 9*x3^2);
                    phi5 = (9/64)*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 - x2)*(1 - x3);
                    phi6 = (9/64)*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 + x2)*(1 - x3);
                    phi7 = (9/64)*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 - x2)*(1 - x3);
                    phi8 = (9/64)*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 + x2)*(1 - x3);
                    phi9 = (1/64)*(1 + x1)*(1 - x2)*(1 - x3)*(-19 + 9*x1^2 + 9*x2^2 + 9*x3^2);
                    phi10 = (9/64)*(1 + x1)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 - x3);
                    phi11 = (9/64)*(1 + x1)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 - x3);
                    phi12 = (1/64)*(1 + x1)*(1 + x2)*(1 - x3)*(-19 + 9*x1^2 + 9*x2^2 + 9*x3^2);
                    phi13 = (9/64)*(1 - x1)*(1 - x2)*(1 - 3*x3 - x3^2 + 3*x3^3);
                    phi14 = (9/64)*(1 - x1)*(1 + x2)*(1 - 3*x3 - x3^2 + 3*x3^3);
                    phi15 = (9/64)*(1 + x1)*(1 - x2)*(1 - 3*x3 - x3^2 + 3*x3^3);
                    phi16 = (9/64)*(1 + x1)*(1 + x2)*(1 - 3*x3 - x3^2 + 3*x3^3);
                    phi17 = (9/64)*(1 - x1)*(1 - x2)*(1 + 3*x3 - x3^2 - 3*x3^3);
                    phi18 = (9/64)*(1 - x1)*(1 + x2)*(1 + 3*x3 - x3^2 - 3*x3^3);
                    phi19 = (9/64)*(1 + x1)*(1 - x2)*(1 + 3*x3 - x3^2 - 3*x3^3);
                    phi20 = (9/64)*(1 + x1)*(1 + x2)*(1 + 3*x3 - x3^2 - 3*x3^3);
                    phi21 = (1/64)*(1 - x1)*(1 - x2)*(1 + x3)*(-19 + 9*x1^2 + 9*x2^2 + 9*x3^2);
                    phi22 = (9/64)*(1 - x1)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 + x3);
                    phi23 = (9/64)*(1 - x1)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 + x3);
                    phi24 = (1/64)*(1 - x1)*(1 + x2)*(1 + x3)*(-19 + 9*x1^2 + 9*x2^2 + 9*x3^2);
                    phi25 = (9/64)*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 - x2)*(1 + x3);
                    phi26 = (9/64)*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 + x2)*(1 + x3);
                    phi27 = (9/64)*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 - x2)*(1 + x3);
                    phi28 = (9/64)*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 + x2)*(1 + x3);
                    phi29 = (1/64)*(1 + x1)*(1 - x2)*(1 + x3)*(-19 + 9*x1^2 + 9*x2^2 + 9*x3^2);
                    phi30 = (9/64)*(1 + x1)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 + x3);
                    phi31 = (9/64)*(1 + x1)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 + x3);
                    phi32 = (1/64)*(1 + x1)*(1 + x2)*(1 + x3)*(-19 + 9*x1^2 + 9*x2^2 + 9*x3^2);                   
                    v = [phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8,... 
                        phi9, phi10, phi11, phi12, phi13, phi14, phi15, phi16,... 
                        phi17, phi18, phi19, phi20, phi21, phi22, phi23, phi24,...
                        phi25, phi26, phi27, phi28, phi29, phi30, phi31, phi32]; 
                end
                
                % Find the indexes of the nodes on the coarse grid for the element (columns of P)
                col = DofsC(elemC,3:3:(3*m));
                
                % Construct the three lines of P corresponding to the degrees of freedom of the node on the fine grid 
                for i = 2:-1:0
                    iP(ind:ind+m-1) = 3*nodef-i;
                    jP(ind:ind+m-1) = col-i;
                    sP(ind:ind+m-1) = v;
                    ind = ind + m;
                end
               
                % Free degrees of freedom for the nodes on the coarse grid 
                if(sum(v==1)==1) % Verify if the node is in the coarse grid 
                    nodec = nodec + 1; % Index of the node on the coarse grid 
                    freedofsc(3*nodec-2) = 3*nodec-2; 
                    freedofsc(3*nodec-1) = 3*nodec-1;
                    freedofsc(3*nodec) = 3*nodec;
                    aux(3*nodec-2) = 3*nodef-2; 
                    aux(3*nodec-1) = 3*nodef-1;
                    aux(3*nodec) = 3*nodef;
                end             
            end          
        end
    end
end

P = sparse(iP(1:ind-1),jP(1:ind-1),sP(1:ind-1));
freedofsc = freedofsc(ismember(aux,freedofsf)); 
P = P(freedofsf,freedofsc); 

end