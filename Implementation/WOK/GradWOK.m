function [Fgrad] = GradWOK(str,u,p,emin,k,x,rmin,w,opfilter,elem,mr,strDens,strDsgn)
% GradWOK calculates the gradient vector of the objective function of the topology optimization problem, without constructing the global stiffness matrix. %
% INPUT: str - structure with the problem data. 
%        u - nodal displacements vector. 
%        p - penalty parameter for the SIMP model. 
%        emin - Young's modulus of the void material
%        k - finite element stiffness matrix. 
%        x - element densities vector (filtered).  
%        rmin - filter radius.
%        w - vector with the sum of the weight factors for each finite element. 
%        opfilter - filter option (0 = no filter, 1 = mean density filter).
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
%        strDens - structure with density elements grid data for multiresolution.
%        strDsgn - structure with design variable elements grid data for multiresolution.
% OUTPUT: Fgrad - gradient vector of the objective function. 
% ---------- % 

if(mr.op)
    Fgrad0 = zeros(strDens.nelem,1);
    nely = strDens.nely; 
    nelxy = strDens.nelx*strDens.nely;   
else 
    Fgrad0 = zeros(str.nelem,1); 
end

elemDeg = elem.deg;

if(elem.type == 1 || elemDeg == 1) % Lagrange elements
    nxny = str.nx*str.ny;
    nxnymnely = nxny-str.ny+1;
    j = 1; % 'j' is the first node of the element 'i'

    aux = zeros(elemDeg+1,elemDeg+1);
    for i2 = 1:(elemDeg+1)
        aux(i2,1) = 1+(i2-1)*nxny;
    end
    for j2 = 2:(elemDeg+1)
        aux(:,j2) = aux(:,1)+(j2-1)*str.ny;
    end
    nodes = zeros((elemDeg+1)^2,(elemDeg+1));
    nodes(:,1) = aux(:);
    for k2 = 2:(elemDeg+1)
        nodes(:,k2) = nodes(:,1)+(k2-1);
    end
    nodes = nodes(:);
    dofs = [3*nodes-2, 3*nodes-1, 3*nodes]; 
    dofs = sort(dofs(:));   

    for el = 1:str.nelem
        uel = u(dofs);
        
        if(mr.op) % Multiresolution
            vb = zeros(mr.n^2, 1); 
            ind = 1; 
            for m1 = 0:(mr.n-1) 
                for m2 = 0:(mr.n-1)
                    vb(ind) = strDens.fde(el) + m1*nelxy + m2*nely; 
                    ind = ind + 1; 
                end
            end     
            v = zeros(mr.n^2, mr.n);
            for m3 = 0:(mr.n-1)
                v(:,m3+1) = vb+m3; 
            end
            v = v(:);         
            v = sort(v);
            
            if (~mr.interp) % Using the original displacements to compute the gradient. 
                for i2 = 1:length(v)  
                    kel = k{i2}; % stiffness integrand for the density element
                    Fgrad0(v(i2)) = -p*(x(v(i2))^(p-1))*(str.E-emin)*(uel'*kel*uel);             
                end
            else % Using interpolated displacements to compute the gradient.
                ud = mr.Idisp*uel; % ud contains the interpolated displacements of all density elements inside the finite element i
                for i2 = 1:length(v) 
                    ueli = ud(mr.Ldofs(i2,:)); % Interpolated displacements of a single density element         
                    kel = k{i2}; % stiffness integrand for the density element
                    Fgrad0(v(i2)) = -p*(x(v(i2))^(p-1))*(str.E-emin)*(ueli'*kel*ueli);             
                end
            end
        else
            Fgrad0(el) = -p*(x(el)^(p-1))*(str.E-emin)*(uel'*k*uel);   
        end
        
        j = j+elemDeg;
        if (mod(j,str.ny) == 0)
            j = j+1+(elemDeg-1)*str.ny;
            if (mod(j,nxny) == nxnymnely)
                j = j+str.ny+(elemDeg-1)*nxny;
                dofs = dofs+(3*(elemDeg+1+elemDeg*str.ny+(elemDeg-1)*nxny));
            else
                dofs = dofs+(3*(elemDeg+1+(elemDeg-1)*str.ny));
            end
        else
            dofs = dofs+(3*elemDeg);
        end
    end
    
elseif(elem.type == 2) % Serendipity elements
    ny1 = str.nely*elemDeg + 1; % number of nodes in the y-direction on layers containing vertices of the elements
    nx2 = str.nelx + 1; % number of nodes in the x direction on layers not containing vertices of the elements 
    ny2 = str.nely + 1; % number of nodes in the y direction on layers not containing vertices of the elements 
    n3 = elemDeg - 1; % number of nodes in the interior of each edge (without counting vertices) of the element
    nxy1 = ny1*(str.nelx+1) + ny2*n3*(str.nelx); % number of nodes in an xy layer containing vertices of the element
    nxy2 = nx2*ny2; % number of nodes in an xy layer not containing vertices of the element

    j = 1; % global index of the first node of the element 
    ex = 1; ey = 1; ez = 1; % layers of the element in each direction 
    for el = 1:str.nelem
        % global indexes of the nodes on vertices of the element
        vt1 = [j, j+elemDeg, j+ny1+n3*ny2, j+ny1+n3*ny2+elemDeg]; 
        vt = [vt1, vt1+nxy1+nxy2*n3]; 

        % global indexes of the nodes on edges of the element in each direction
        edx = zeros(1,4*n3);
        edy = zeros(1,4*n3); 
        edz = zeros(1,4*n3); 
        ind = 1;
        for i = 1:n3
            v1 = vt(1) + ((str.nely-ey+1)*elemDeg+ey) + (i-1)*ny2; % from the vertex vt(1) go in the x-direction to the node in the next layer (v1)
            v2 = vt(5) + ((str.nely-ey+1)*elemDeg+ey) + (i-1)*ny2; % from the vertex vt(5) go in the x-direction to the node in the next layer (v2)
            v3 = vt(1) + (nxy1-((ex-1)*(ny1+n3*ny2)) + (ex-1)*ny2 + (ey-1)*(1-elemDeg)) + (i-1)*nxy2; % from the vertex vt(1) go in the z-direction to the node v3
            v4 = v3+ny2; % from the vertex vt(3) go in the z-direction to the node v4 (or from v3 go in the x-direction)
            edx(ind:ind+3) = [v1, v1+1, v2, v2+1]; 
            edy(ind:ind+3) = vt([1,3,5,7]) + i; 
            edz(ind:ind+3) = [v3, v3+1, v4, v4+1];
            ind = ind + 4; 
        end
 
        nodes = [vt,edx,edy,edz]; % indexes of all nodes of the element 
        dofs = [3*nodes-2, 3*nodes-1, 3*nodes]; % degrees of freedom of the element
        dofs = sort(dofs(:));
        
        uel = u(dofs);
        
        if(mr.op) % Multiresolution
            vb = zeros(mr.n^2, 1); 
            ind = 1; 
            for m1 = 0:(mr.n-1) 
                for m2 = 0:(mr.n-1)
                    vb(ind) = strDens.fde(el) + m1*nelxy + m2*nely; 
                    ind = ind + 1; 
                end
            end     
            v = zeros(mr.n^2, mr.n);
            for m3 = 0:(mr.n-1)
                v(:,m3+1) = vb+m3; 
            end
            v = v(:);         
            v = sort(v);
            
            if (~mr.interp) % Using the original displacements to compute the gradient. 
                for i2 = 1:length(v)  
                    kel = k{i2}; % stiffness integrand for the density element
                    Fgrad0(v(i2)) = -p*(x(v(i2))^(p-1))*(str.E-emin)*(uel'*kel*uel);             
                end
            else % Using interpolated displacements to compute the gradient.
                ud = mr.Idisp*uel; % ud contains the interpolated displacements of all density elements inside the finite element i
                for i2 = 1:length(v) 
                    ueli = ud(mr.Ldofs(i2,:)); % Interpolated displacements of a single density element         
                    kel = k{i2}; % stiffness integrand for the density element
                    Fgrad0(v(i2)) = -p*(x(v(i2))^(p-1))*(str.E-emin)*(ueli'*kel*ueli);             
                end
            end
        else
            Fgrad0(el) = -p*(x(el)^(p-1))*(str.E-emin)*(uel'*k*uel); 
        end

        % Update the layers in each direction and the index of the first node for the next element
        if(ey == str.nely)
            ey = 1; 
            if(ex == str.nelx)
                ex = 1; 
                ez = ez + 1;
                j = j + elemDeg + 1 + n3*ny2 + ny1 + n3*nxy2; 
            else
                ex = ex + 1; 
                j = j + elemDeg + 1 + n3*ny2;
            end
        else
            ey = ey + 1; 
            j = j + elemDeg;
        end
    end
end

Fgrad = ApplyFilterWOK(Fgrad0,str,rmin,w,opfilter,true,mr,strDens,strDsgn);

end