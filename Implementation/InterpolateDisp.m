function [Idisp,Ldofs,kd] = InterpolateDisp(mr,elem,strDens)
% InterpolateDisp approximates the nodal displacements for a density element, using the finite element displacements, for multiresolution. %
% INPUT: mr - structure that contains parameters used by the multiresolution method.
%        elem - structure with element characteristics (deg - polynomial degree, type - element type (1 = Lagrange, 2 = serendipity).
%        strDens - structure with density mesh data for multiresolution.
% OUTPUT: Idisp - interpolation matrix to obtain the nodal displacements of the density elements inside a finite element.
%         Ldofs - matrix with the local degrees of freedom of each density element inside a finite element. 
%         kd - stiffness matrix of the density element. 
% ---------- %

elemDeg = elem.deg;

% m = total number of nodes in each element
if(elem.type == 1) % Lagrange element
    m = (elemDeg+1)^3;
elseif(elem.type == 2) % Serendipity element
    m = 8 + 12*(elemDeg-1);
end

% n = total dofs of the density elements inside each finite element
if(elem.type == 1 || elemDeg == 1) % Lagrange element
    n = 3*((mr.n*elemDeg+1)^3);
elseif(elem.type == 2) % Serendipity element
    nxy1 = (mr.n*elemDeg+1)*(mr.n+1) + (elemDeg-1)*(mr.n+1)*(mr.n);
    n = 3*(nxy1*(mr.n+1) + (elemDeg-1)*((mr.n+1)^2)*(mr.n));
end

iIdisp = zeros(n,1);
jIdisp = zeros(n,1);
sIdisp = zeros(n,1);

row = 1;
col = 1:3:(3*m);
ind = 1;
indstep = m - 1;
cordstep = 2/(mr.n*elemDeg);
nt = [1,1,1]; % local layers of the node inside the element
for x3 = -1:cordstep:1
    for x1 = -1:cordstep:1
        for x2 = -1:cordstep:1
            if(elem.type == 1 || elemDeg == 1 || (elem.type == 2 && sum(mod(nt,elemDeg)==1) >= 2)) % verify if the point nt is a node for the serendipity element 
                if(elemDeg == 1) % Linear element (8 shape functions)
                    phi1 = (1/8)*(1-x1)*(1-x2)*(1-x3);
                    phi2 = (1/8)*(1-x1)*(1+x2)*(1-x3);
                    phi3 = (1/8)*(1+x1)*(1-x2)*(1-x3);
                    phi4 = (1/8)*(1+x1)*(1+x2)*(1-x3);
                    phi5 = (1/8)*(1-x1)*(1-x2)*(1+x3);
                    phi6 = (1/8)*(1-x1)*(1+x2)*(1+x3);
                    phi7 = (1/8)*(1+x1)*(1-x2)*(1+x3);
                    phi8 = (1/8)*(1+x1)*(1+x2)*(1+x3);
                    v = [phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8]; % Shape functions
                elseif(elemDeg == 2)
                    if(elem.type == 1) % Lagrange quadratic element (27 shape functions)
                        phi1 = (1/8)*(x1^2 - x1)*(x2^2 - x2)*(x3^2 - x3);
                        phi2 = (1/4)*(x1^2 - x1)*(1 - x2^2)*(x3^2 - x3);
                        phi3 = (1/8)*(x1^2 - x1)*(x2^2 + x2)*(x3^2 - x3);
                        phi4 = (1/4)*(1 - x1^2)*(x2^2 - x2)*(x3^2 - x3);
                        phi5 = (1/2)*(1 - x1^2)*(1 - x2^2)*(x3^2 - x3);
                        phi6 = (1/4)*(1 - x1^2)*(x2^2 + x2)*(x3^2 - x3);
                        phi7 = (1/8)*(x1^2 + x1)*(x2^2 - x2)*(x3^2 - x3);
                        phi8 = (1/4)*(x1^2 + x1)*(1 - x2^2)*(x3^2 - x3);
                        phi9 = (1/8)*(x1^2 + x1)*(x2^2 + x2)*(x3^2 - x3);
                        phi10 = (1/4)*(x1^2 - x1)*(x2^2 - x2)*(1 - x3^2);
                        phi11 = (1/2)*(x1^2 - x1)*(1 - x2^2)*(1 - x3^2);
                        phi12 = (1/4)*(x1^2 - x1)*(x2^2 + x2)*(1 - x3^2);
                        phi13 = (1/2)*(1 - x1^2)*(x2^2 - x2)*(1 - x3^2);
                        phi14 = (1 - x1^2)*(1 - x2^2)*(1 - x3^2);
                        phi15 = (1/2)*(1 - x1^2)*(x2^2 + x2)*(1 - x3^2);
                        phi16 = (1/4)*(x1^2 + x1)*(x2^2 - x2)*(1 - x3^2);
                        phi17 = (1/2)*(x1^2 + x1)*(1 - x2^2)*(1 - x3^2);
                        phi18 = (1/4)*(x1^2 + x1)*(x2^2 + x2)*(1 - x3^2);
                        phi19 = (1/8)*(x1^2 - x1)*(x2^2 - x2)*(x3^2 + x3);
                        phi20 = (1/4)*(x1^2 - x1)*(1 - x2^2)*(x3^2 + x3);
                        phi21 = (1/8)*(x1^2 - x1)*(x2^2 + x2)*(x3^2 + x3);
                        phi22 = (1/4)*(1 - x1^2)*(x2^2 - x2)*(x3^2 + x3);
                        phi23 = (1/2)*(1 - x1^2)*(1 - x2^2)*(x3^2 + x3);
                        phi24 = (1/4)*(1 - x1^2)*(x2^2 + x2)*(x3^2 + x3);
                        phi25 = (1/8)*(x1^2 + x1)*(x2^2 - x2)*(x3^2 + x3);
                        phi26 = (1/4)*(x1^2 + x1)*(1 - x2^2)*(x3^2 + x3);
                        phi27 = (1/8)*(x1^2 + x1)*(x2^2 + x2)*(x3^2 + x3);
                        v = [phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8, phi9,...
                            phi10, phi11, phi12, phi13, phi14, phi15, phi16, phi17, phi18, ...
                            phi19, phi20, phi21, phi22, phi23, phi24, phi25, phi26, phi27]; 
                    elseif(elem.type == 2) % Serendipity quadratic element (20 shape functions)
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
                    end
                elseif(elemDeg == 3)
                    if(elem.type == 1) % Lagrange cubic element (64 shape functions)                      
                        phi1 = 1/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi2 = 9/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi3 = 9/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi4 = 1/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi5 = 9/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi6 = 81/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi7 = 81/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi8 = 9/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3); 
                        phi9 = 9/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi10 = 81/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3); 
                        phi11 = 81/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi12 = 9/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi13 = 1/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi14 = 9/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi15 = 9/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi16 = 1/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(-1 + x3 + 9*x3^2 - 9*x3^3);
                        phi17 = 9/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi18 = 81/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi19 = 81/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi20 = 9/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi21 = 81/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi22 = 729/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi23 = 729/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi24 = 81/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi25 = 81/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi26 = 729/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3); 
                        phi27 = 729/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi28 = 81/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi29 = 9/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi30 = 81/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi31 = 81/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi32 = 9/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(1 - 3*x3 - x3^2 + 3*x3^3);
                        phi33 = 9/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi34 = 81/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi35 = 81/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3); 
                        phi36 = 9/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi37 = 81/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3); 
                        phi38 = 729/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi39 = 729/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi40 = 81/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi41 = 81/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi42 = 729/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi43 = 729/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi44 = 81/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi45 = 9/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi46 = 81/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi47 = 81/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi48 = 9/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(1 + 3*x3 - x3^2 - 3*x3^3);
                        phi49 = 1/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi50 = 9/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi51 = 9/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi52 = 1/4096*(-1 + x1 + 9*x1^2 - 9*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi53 = 9/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi54 = 81/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi55 = 81/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi56 = 9/4096*(1 - 3*x1 - x1^2 + 3*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi57 = 9/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi58 = 81/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi59 = 81/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi60 = 9/4096*(1 + 3*x1 - x1^2 - 3*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi61 = 1/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(-1 + x2 + 9*x2^2 - 9*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi62 = 9/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(1 - 3*x2 - x2^2 + 3*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi63 = 9/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(1 + 3*x2 - x2^2 - 3*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        phi64 = 1/4096*(-1 - x1 + 9*x1^2 + 9*x1^3)*(-1 - x2 + 9*x2^2 + 9*x2^3)*(-1 - x3 + 9*x3^2 + 9*x3^3);
                        v = [phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8,... 
                            phi9, phi10, phi11, phi12, phi13, phi14, phi15, phi16,... 
                            phi17, phi18, phi19, phi20, phi21, phi22, phi23, phi24,...
                            phi25, phi26, phi27, phi28, phi29, phi30, phi31, phi32,...
                            phi33, phi34, phi35, phi36, phi37, phi38, phi39, phi40,...
                            phi41, phi42, phi43, phi44, phi45, phi46, phi47, phi48,...
                            phi49, phi50, phi51, phi52, phi53, phi54, phi55, phi56,...
                            phi57, phi58, phi59, phi60, phi61, phi62, phi63, phi64]; 
                    elseif(elem.type == 2) % Serendipity cubic element (32 shape functions)
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
                end

                for i = 0:2
                    iIdisp(ind:ind+indstep) = row;
                    jIdisp(ind:ind+indstep) = col + i;
                    sIdisp(ind:ind+indstep) = v;
                    ind = ind+indstep+1;
                    row = row+1;
                end
            end
            nt(2) = nt(2) + 1;
        end
        nt(1) = nt(1) + 1;
        nt(2) = 1;
    end
    nt(3) = nt(3) + 1;
    nt(1) = 1;
end

Idisp = sparse(iIdisp,jIdisp,sIdisp,n,3*m); % Interpolation matrix

% Obtaining the local degrees of freedom of each density element inside a finite element
if(elem.type == 1 || elemDeg == 1)
    strL.nx = mr.n*elemDeg + 1;
    strL.ny = strL.nx;
    strL.nelem = (mr.n)^3;
elseif(elem.type == 2) 
    strL.nelx = mr.n;
    strL.nely = mr.n;
    strL.nelz = mr.n;
    strL.nelem = (mr.n)^3;
end
Ldofs = DegreesOfFreedom(strL,elem); 

% Element stiffness matrix of the density element
kd = ElemStiffMatrixGQ(strDens,elem.deg+1,elem); 

end