function y = KelemProd(str,x,kv,dens,p,emin,elem,mr)
% KelemProd calculates the product of the global stiffness matrix K by a vector, without constructing the matrix K. %
% INPUT: str - structure with the problem data. 
%        x - function argument (vector that will be multiplied by K).
%        kv - element stiffness matrix, converted into a column vector.
%        dens - element densities vector (filtered). 
%        p - penalty parameter for the SIMP model. 
%        emin - Young's modulus of the void material
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: y - product K*x.
% ---------- %

y = zeros(str.nnodes*3,1);
xc = zeros(str.nnodes*3,1);
xc(str.freedofs) = x;

elemDeg = elem.deg;

% m = total number of nodes in each element
if(elem.type == 1) % Lagrange element
    m = (elemDeg+1)^3;
elseif(elem.type == 2) % Serendipity element (with degree 1, 2 or 3)
    m = 8 + 12*(elemDeg-1);
end

if(mr.op)
    jmr = 1; 
    nelx = str.nelx*(mr.n);
    nely = str.nely*(mr.n);
    nelxy = nelx*nely;
else
    kel = reshape(kv,[3*m,3*m]);
end

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
        if(mr.op) % Multiresolution
            vb = zeros(mr.n^2, 1); 
            ind = 1; 
            for m1 = 0:(mr.n-1) 
                for m2 = 0:(mr.n-1)
                    vb(ind) = jmr + m1*nelxy + m2*nely; 
                    ind = ind + 1; 
                end
            end         
            v = zeros(mr.n^2, mr.n);
            for m3 = 0:(mr.n-1)
                v(:,m3+1) = vb+m3; 
            end
            v = v(:); 
            v = sort(v);

            kel = kv*(emin+(str.E-emin)*(dens(v)'.^p))';
            y(dofs) = y(dofs)+reshape(kel,[3*m,3*m])*xc(dofs);
            
            jmr = jmr+mr.n;
            if (mod(jmr-1,nely) == 0)
                jmr = jmr+(mr.n-1)*nely;
                if (mod(jmr,nelxy) == 1)
                    jmr = jmr+(mr.n-1)*nelxy;
                end
            end
        else
            y(dofs) = y(dofs)+kel*((emin+(str.E-emin)*(dens(el)^p))*xc(dofs));
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
        dofs = sort(dofs); 
        
        if(mr.op) % Multiresolution
            vb = zeros(mr.n^2, 1); 
            ind = 1; 
            for m1 = 0:(mr.n-1) 
                for m2 = 0:(mr.n-1)
                    vb(ind) = jmr + m1*nelxy + m2*nely; 
                    ind = ind + 1; 
                end
            end         
            v = zeros(mr.n^2, mr.n);
            for m3 = 0:(mr.n-1)
                v(:,m3+1) = vb+m3; 
            end
            v = v(:); 
            v = sort(v);

            kel = kv*(emin+(str.E-emin)*(dens(v)'.^p))';
            y(dofs) = y(dofs)+reshape(kel,[3*m,3*m])*xc(dofs);
            
            jmr = jmr+mr.n;
            if (mod(jmr-1,nely) == 0)
                jmr = jmr+(mr.n-1)*nely;
                if (mod(jmr,nelxy) == 1)
                    jmr = jmr+(mr.n-1)*nelxy;
                end
            end
        else
            y(dofs) = y(dofs)+kel*((emin+(str.E-emin)*(dens(el)^p))*xc(dofs));
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

y = y(str.freedofs);

end