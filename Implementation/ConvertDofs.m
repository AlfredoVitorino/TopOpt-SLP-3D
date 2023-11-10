function [str] = ConvertDofs(str,elem)
% ConvertDofs adapts the number of nodes, freedofs and the vector of nodal loads according to the polynomial degree of the shape functions. %
% INPUT: str - structure with the problem data. 
%        elem - structure with element characteristics (deg - polynomial degree, type - element type (1 = Lagrange, 2 = serendipity).
% OUTPUT: str - structure with the problem data updated.
% ---------- %

elemDeg = elem.deg;

if(elem.type == 1 || elemDeg == 1) % Lagrange element  
    str.nx = str.nelx*elemDeg+1; % new number of nodes in the x direction
    str.ny = str.nely*elemDeg+1; % new number of nodes in the y direction
    str.nz = str.nelz*elemDeg+1; % new number of nodes in the z direction
    str.nnodes = str.nx*str.ny*str.nz; % new total number of nodes
    elnew = str.el/elemDeg; % distance between two nodes in the x direction 
    ehnew = str.eh/elemDeg; % distance between two nodes in the y direction 
    ewnew = str.ew/elemDeg; % distance between two nodes in the z direction 
    
    % Converting supports %
    j = 1;
    for i = 1:size(str.supp,1)
        % Initial and final coordinates of the support
        x0 = str.supp(i,1);
        xf = str.supp(i,2);
        y0 = str.supp(i,3);
        yf = str.supp(i,4);
        z0 = str.supp(i,5); 
        zf = str.supp(i,6);
        % Indicates if the support prevents the movement in each direction
        ix = str.supp(i,7); 
        iy = str.supp(i,8); 
        iz = str.supp(i,9);
        
        nx0 = 1 + floor(x0/elnew) + (mod(x0,elnew) > elnew/2); % Initial nodes layer of the support in the x direction 
        nxf = 1 + floor(xf/elnew) + (mod(xf,elnew) > elnew/2); % Final nodes layer of the support in the x direction
        ny0 = 1 + floor(y0/ehnew) + (mod(y0,ehnew) > ehnew/2); % Initial nodes layer of the support in the y direction
        nyf = 1 + floor(yf/ehnew) + (mod(yf,ehnew) > ehnew/2); % Final nodes layer of the support in the y direction
        nz0 = 1 + floor(z0/ewnew) + (mod(z0,ewnew) > ewnew/2); % Initial nodes layer of the support in the z direction
        nzf = 1 + floor(zf/ewnew) + (mod(zf,ewnew) > ewnew/2); % Final nodes layer of the support in the z direction
        
        node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Initial node of the support
        for nz = nz0:nzf
            for nx = nx0:nxf
                for ny = ny0:nyf
                    supp.nodes(j) = node;
                    supp.ix(j) = ix;
                    supp.iy(j) = iy;
                    supp.iz(j) = iz;
                    node = node + 1;
                    j = j + 1;
                end
                node = node - 1 + (str.ny - nyf) + ny0;
            end
            node = node + str.ny*(str.nx - nxf + nx0 - 1);
        end
    end
    s = 3*supp.nodes;
    ind1 = [s-2; s-1; s];
    ind1 = ind1(:);
    ind2 = [supp.ix; supp.iy; supp.iz]; 
    ind2 = ind2(:);
    supp.ind = ind1(ind2~=0); % Global indexes of the nodes with supports preventing the movement
    str.freedofs = 1:3*str.nnodes; 
    str.freedofs(supp.ind) = []; % Indexes of nonzero degrees of freedom. 
    % ---------- %
    
    % Converting loads %
    j = 1;
    for i = 1:size(str.loads,1)
        % Initial and final coordinates of the load application
        x0 = str.loads(i,1);
        xf = str.loads(i,2);
        y0 = str.loads(i,3);
        yf = str.loads(i,4);
        z0 = str.loads(i,5);
        zf = str.loads(i,6);
        % Components of load
        Fx = str.loads(i,7); 
        Fy = str.loads(i,8); 
        Fz = str.loads(i,9);
        
        nx0 = 1 + floor(x0/elnew) + (mod(x0,elnew) > elnew/2); % Initial nodes layer of load application in the x direction 
        nxf = 1 + floor(xf/elnew) + (mod(xf,elnew) > elnew/2); % Final nodes layer of load application in the x direction 
        ny0 = 1 + floor(y0/ehnew) + (mod(y0,ehnew) > ehnew/2); % Initial nodes layer of load application in the y direction
        nyf = 1 + floor(yf/ehnew) + (mod(yf,ehnew) > ehnew/2); % Final nodes layer of load application in the y direction
        nz0 = 1 + floor(z0/ewnew) + (mod(z0,ewnew) > ewnew/2); % Initial nodes layer of load application in the z direction
        nzf = 1 + floor(zf/ewnew) + (mod(zf,ewnew) > ewnew/2); % Final nodes layer of load application in the z direction

        node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Initial node where the load is applied
        
        % Concentrated Load %
        if(x0 == xf && y0 == yf && z0 == zf) 
            loads.nodes(j) = node;
            loads.Fx(j) = Fx;
            loads.Fy(j) = Fy;
            loads.Fz(j) = Fz;
            j = j+1; 
            continue;
        end
        
        % Load distributed only in the x direction %
        if(x0 ~= xf && y0 == yf && z0 == zf) 
            if(nx0 == nxf)
                loads.nodes(j) = node; 
                loads.Fx(j) = (xf-x0)*Fx;
                loads.Fy(j) = (xf-x0)*Fy;
                loads.Fz(j) = (xf-x0)*Fz;
                j = j+1;
            else
                loads.nodes(j) = node; 
                loads.Fx(j) = ((nx0-1/2)*elnew - x0)*Fx;
                loads.Fy(j) = ((nx0-1/2)*elnew - x0)*Fy;
                loads.Fz(j) = ((nx0-1/2)*elnew - x0)*Fz;
                j = j+1;
                node = node + str.ny; 
                for i2 = (nx0+1):(nxf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = elnew*Fx;
                    loads.Fy(j) = elnew*Fy;
                    loads.Fz(j) = elnew*Fz;
                    j = j+1; 
                    node = node + str.ny;
                end
                loads.nodes(j) = node; 
                loads.Fx(j) = (xf - (nxf-3/2)*elnew)*Fx;
                loads.Fy(j) = (xf - (nxf-3/2)*elnew)*Fy;
                loads.Fz(j) = (xf - (nxf-3/2)*elnew)*Fz;
                j = j+1;
            end
            continue;
        end
        
        % Load distributed only in the y direction %
        if(x0 == xf && y0 ~= yf && z0 == zf)
            if(ny0 == nyf)
                loads.nodes(j) = node; 
                loads.Fx(j) = (yf-y0)*Fx;
                loads.Fy(j) = (yf-y0)*Fy;
                loads.Fz(j) = (yf-y0)*Fz;
                j = j+1;
            else
                loads.nodes(j) = node; 
                loads.Fx(j) = ((ny0-1/2)*ehnew - y0)*Fx;
                loads.Fy(j) = ((ny0-1/2)*ehnew - y0)*Fy;
                loads.Fz(j) = ((ny0-1/2)*ehnew - y0)*Fz;
                j = j+1;
                node = node + 1; 
                for i2 = (ny0+1):(nyf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = ehnew*Fx;
                    loads.Fy(j) = ehnew*Fy;
                    loads.Fz(j) = ehnew*Fz;
                    j = j+1; 
                    node = node + 1;
                end
                loads.nodes(j) = node; 
                loads.Fx(j) = (yf - (nyf-3/2)*ehnew)*Fx;
                loads.Fy(j) = (yf - (nyf-3/2)*ehnew)*Fy;
                loads.Fz(j) = (yf - (nyf-3/2)*ehnew)*Fz;
                j = j+1;
            end
            continue;
        end
        
        % Load distributed only in the z direction %
        if(x0 == xf && y0 == yf && z0 ~= zf)
            if(nz0 == nzf)
                loads.nodes(j) = node; 
                loads.Fx(j) = (zf-z0)*Fx;
                loads.Fy(j) = (zf-z0)*Fy;
                loads.Fz(j) = (zf-z0)*Fz;
                j = j+1;
            else
                loads.nodes(j) = node; 
                loads.Fx(j) = ((nz0-1/2)*ewnew - z0)*Fx;
                loads.Fy(j) = ((nz0-1/2)*ewnew - z0)*Fy;
                loads.Fz(j) = ((nz0-1/2)*ewnew - z0)*Fz;
                j = j+1;
                node = node + str.nx*str.ny; 
                for i2 = (nz0+1):(nzf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = ewnew*Fx;
                    loads.Fy(j) = ewnew*Fy;
                    loads.Fz(j) = ewnew*Fz;
                    j = j+1; 
                    node = node + str.nx*str.ny;
                end
                loads.nodes(j) = node; 
                loads.Fx(j) = (zf - (nzf-3/2)*ewnew)*Fx;
                loads.Fy(j) = (zf - (nzf-3/2)*ewnew)*Fy;
                loads.Fz(j) = (zf - (nzf-3/2)*ewnew)*Fz;
                j = j+1;
            end
            continue;
        end
        
        % Load distributed in an area parallel to the xy plane %
        if(x0 ~= xf && y0 ~= yf && z0 == zf)
            if(nx0 == nxf)
                d = xf-x0;
                if(ny0 == nyf)
                    loads.nodes(j) = node;
                    loads.Fx(j) = d*(yf-y0)*Fx; 
                    loads.Fy(j) = d*(yf-y0)*Fy;
                    loads.Fz(j) = d*(yf-y0)*Fz;
                    j = j+1; 
                else
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    node = node + 1; 
                    for i2 = (ny0+1):(nyf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                        node = node + 1;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyf-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyf-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyf-3/2)*ehnew)*Fz;
                    j = j+1;
                end 
            else
                if(ny0 == nyf)
                    d = yf-y0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    node = node + str.ny; 
                    for i2 = (nx0+1):(nxf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                        node = node + str.ny;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxf-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxf-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxf-3/2)*elnew)*Fz;
                    j = j+1;
                else 
                    d = (nx0-1/2)*elnew - x0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    node = node + 1; 
                    for i2 = (ny0+1):(nyf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                        node = node + 1;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyf-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyf-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyf-3/2)*ehnew)*Fz;
                    j = j+1;
                    node = node + (str.ny - nyf) + ny0;
                    
                    d = elnew;
                    for j2 = (nx0+1):(nxf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*((ny0-1/2)*ehnew - y0)*Fx;
                        loads.Fy(j) = d*((ny0-1/2)*ehnew - y0)*Fy;
                        loads.Fz(j) = d*((ny0-1/2)*ehnew - y0)*Fz;
                        j = j+1;
                        node = node + 1; 
                        for i2 = (ny0+1):(nyf-1)
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ehnew*Fx;
                            loads.Fy(j) = d*ehnew*Fy;
                            loads.Fz(j) = d*ehnew*Fz;
                            j = j+1; 
                            node = node + 1;
                        end
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*(yf - (nyf-3/2)*ehnew)*Fx;
                        loads.Fy(j) = d*(yf - (nyf-3/2)*ehnew)*Fy;
                        loads.Fz(j) = d*(yf - (nyf-3/2)*ehnew)*Fz;
                        j = j+1;
                        node = node + (str.ny - nyf) + ny0;
                    end
                    
                    d = xf - (nxf-3/2)*elnew;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    node = node + 1; 
                    for i2 = (ny0+1):(nyf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                        node = node + 1;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyf-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyf-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyf-3/2)*ehnew)*Fz;
                    j = j+1;
                end 
            end
            continue;
        end
        
        % Load distributed in an area parallel to the xz plane %
        if(x0 ~= xf && y0 == yf && z0 ~= zf)
            if(nx0 == nxf)
                d = xf-x0;
                if(nz0 == nzf)
                    loads.nodes(j) = node;
                    loads.Fx(j) = d*(zf-z0)*Fx; 
                    loads.Fy(j) = d*(zf-z0)*Fy;
                    loads.Fz(j) = d*(zf-z0)*Fz;
                    j = j+1; 
                else
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nz0-1/2)*ewnew - z0)*Fx;
                    loads.Fy(j) = d*((nz0-1/2)*ewnew - z0)*Fy;
                    loads.Fz(j) = d*((nz0-1/2)*ewnew - z0)*Fz;
                    j = j+1;
                    node = node + str.nx*str.ny; 
                    for i2 = (nz0+1):(nzf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ewnew*Fx;
                        loads.Fy(j) = d*ewnew*Fy;
                        loads.Fz(j) = d*ewnew*Fz;
                        j = j+1; 
                        node = node + str.nx*str.ny;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(zf - (nzf-3/2)*ewnew)*Fx;
                    loads.Fy(j) = d*(zf - (nzf-3/2)*ewnew)*Fy;
                    loads.Fz(j) = d*(zf - (nzf-3/2)*ewnew)*Fz;
                    j = j+1;
                end 
            else
                if(nz0 == nzf)
                    d = zf-z0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    node = node + str.ny; 
                    for i2 = (nx0+1):(nxf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                        node = node + str.ny;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxf-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxf-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxf-3/2)*elnew)*Fz;
                    j = j+1;
                else 
                    d = (nz0-1/2)*ewnew - z0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    node = node + str.ny; 
                    for i2 = (nx0+1):(nxf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                        node = node + str.ny;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxf-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxf-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxf-3/2)*elnew)*Fz;
                    j = j+1;
                    node = node + str.nx*str.ny - (nxf-nx0)*str.ny;

                    d = ewnew;
                    for j2 = (nz0+1):(nzf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*((nx0-1/2)*elnew - x0)*Fx;
                        loads.Fy(j) = d*((nx0-1/2)*elnew - x0)*Fy;
                        loads.Fz(j) = d*((nx0-1/2)*elnew - x0)*Fz;
                        j = j+1;
                        node = node + str.ny; 
                        for i2 = (nx0+1):(nxf-1)
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                            node = node + str.ny;
                        end
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*(xf - (nxf-3/2)*elnew)*Fx;
                        loads.Fy(j) = d*(xf - (nxf-3/2)*elnew)*Fy;
                        loads.Fz(j) = d*(xf - (nxf-3/2)*elnew)*Fz;
                        j = j+1;
                        node = node + str.nx*str.ny - (nxf-nx0)*str.ny;
                    end

                    d = zf - (nzf-3/2)*ewnew;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    node = node + str.ny; 
                    for i2 = (nx0+1):(nxf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                        node = node + str.ny;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxf-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxf-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxf-3/2)*elnew)*Fz;
                    j = j+1;
                end 
            end
            continue;
        end
        
        % Load distributed in an area parallel to the yz plane %
        if(x0 == xf && y0 ~= yf && z0 ~= zf)
            if(nz0 == nzf)
                d = zf-z0;
                if(ny0 == nyf)
                    loads.nodes(j) = node;
                    loads.Fx(j) = d*(yf-y0)*Fx; 
                    loads.Fy(j) = d*(yf-y0)*Fy;
                    loads.Fz(j) = d*(yf-y0)*Fz;
                    j = j+1; 
                else
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    node = node + 1; 
                    for i2 = (ny0+1):(nyf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                        node = node + 1;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyf-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyf-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyf-3/2)*ehnew)*Fz;
                    j = j+1;
                end 
            else
                if(ny0 == nyf)
                    d = yf-y0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nz0-1/2)*ewnew - z0)*Fx;
                    loads.Fy(j) = d*((nz0-1/2)*ewnew - z0)*Fy;
                    loads.Fz(j) = d*((nz0-1/2)*ewnew - z0)*Fz;
                    j = j+1;
                    node = node + str.nx*str.ny; 
                    for i2 = (nz0+1):(nzf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ewnew*Fx;
                        loads.Fy(j) = d*ewnew*Fy;
                        loads.Fz(j) = d*ewnew*Fz;
                        j = j+1; 
                        node = node + str.nx*str.ny;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(zf - (nzf-3/2)*ewnew)*Fx;
                    loads.Fy(j) = d*(zf - (nzf-3/2)*ewnew)*Fy;
                    loads.Fz(j) = d*(zf - (nzf-3/2)*ewnew)*Fz;
                    j = j+1;
                else 
                    d = (nz0-1/2)*ewnew - z0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    node = node + 1; 
                    for i2 = (ny0+1):(nyf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                        node = node + 1;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyf-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyf-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyf-3/2)*ehnew)*Fz;
                    j = j+1;
                    node = node + str.nx*str.ny - (nyf-ny0);

                    d = ewnew;
                    for j2 = (nz0+1):(nzf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*((ny0-1/2)*ehnew - y0)*Fx;
                        loads.Fy(j) = d*((ny0-1/2)*ehnew - y0)*Fy;
                        loads.Fz(j) = d*((ny0-1/2)*ehnew - y0)*Fz;
                        j = j+1;
                        node = node + 1; 
                        for i2 = (ny0+1):(nyf-1)
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ehnew*Fx;
                            loads.Fy(j) = d*ehnew*Fy;
                            loads.Fz(j) = d*ehnew*Fz;
                            j = j+1; 
                            node = node + 1;
                        end
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*(yf - (nyf-3/2)*ehnew)*Fx;
                        loads.Fy(j) = d*(yf - (nyf-3/2)*ehnew)*Fy;
                        loads.Fz(j) = d*(yf - (nyf-3/2)*ehnew)*Fz;
                        j = j+1;
                        node = node + str.nx*str.ny - (nyf-ny0);
                    end

                    d = zf - (nzf-3/2)*ewnew;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    node = node + 1; 
                    for i2 = (ny0+1):(nyf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                        node = node + 1;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyf-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyf-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyf-3/2)*ehnew)*Fz;
                    j = j+1;
                end 
            end
            continue;
        end  
    end
    
    str.f = zeros(3*str.nnodes, 1); % Vector of nodal loads
    s = 3*loads.nodes; 
    for i = 1:length(s)
        str.f(s(i)-2) = str.f(s(i)-2) + loads.Fx(i); 
        str.f(s(i)-1) = str.f(s(i)-1) + loads.Fy(i);
        str.f(s(i)) = str.f(s(i)) + loads.Fz(i);
    end 
    str.f(supp.ind) = []; % Boundary conditions
    % ---------- %
    
elseif(elem.type == 2) % serendipity element
    ny1 = str.nely*elemDeg + 1; % number of nodes in the y direction on layers containing vertices of the elements
    nx2 = str.nelx + 1; % number of nodes in the x direction on layers not containing vertices of the elements 
    ny2 = str.nely + 1; % number of nodes in the y direction on layers not containing vertices of the elements 
    n3 = elemDeg - 1; % number of nodes in the interior of each edge (without counting vertices) of the element 
    nxy1 = ny1*(str.nelx+1) + ny2*n3*(str.nelx); % number of nodes in an xy layer containing vertices of the element
    nxy2 = nx2*ny2; % number of nodes in an xy layer not containing vertices of the element
    nnodes = nxy1*(str.nelz+1) + nxy2*n3*str.nelz;  % total number of nodes on the grid   
    elnew = str.el/elemDeg; % distance between two nodes in the x direction on layers containing vertices of the elements 
    ehnew = str.eh/elemDeg; % distance between two nodes in the y direction on layers containing vertices of the elements
    ewnew = str.ew/elemDeg; % distance between two nodes in the z direction on layers containing vertices of the elements

    % Converting supports %
    j = 1;
    for i = 1:size(str.supp,1)
        % Initial and final coordinates of the support
        x0 = str.supp(i,1);
        xf = str.supp(i,2);
        y0 = str.supp(i,3);
        yf = str.supp(i,4);
        z0 = str.supp(i,5); 
        zf = str.supp(i,6);
        % Indicates if the support prevents the movement in each direction
        ix = str.supp(i,7); 
        iy = str.supp(i,8); 
        iz = str.supp(i,9);
        
        nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Initial nodes layer (on the linear grid) of the support in the x direction 
        nxf = 1 + floor(xf/str.el) + (mod(xf,str.el) > str.el/2); % Final nodes layer (on the linear grid) of the support in the x direction
        ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Initial nodes layer (on the linear grid) of the support in the y direction
        nyf = 1 + floor(yf/str.eh) + (mod(yf,str.eh) > str.eh/2); % Final nodes layer (on the linear grid) of the support in the y direction
        nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Initial nodes layer (on the linear grid) of the support in the z direction
        nzf = 1 + floor(zf/str.ew) + (mod(zf,str.ew) > str.ew/2); % Final nodes layer (on the linear grid) of the support in the z direction

        for nz = nz0:nzf
            for nx = nx0:nxf
                for ny = ny0:nyf
                    node = (nz-1)*(nxy1+n3*nxy2) + (nx-1)*(ny1+n3*ny2) + (ny-1)*(1+n3) + 1;
                    if(ny == nyf)
                        suppf.nodes(j) = node; 
                        suppf.ix(j) = ix; 
                        suppf.iy(j) = iy;
                        suppf.iz(j) = iz;
                        j = j+1;
                    else % add the nodes in the interior of the edges
                        for j2 = 0:n3
                            suppf.nodes(j) = node+j2; 
                            suppf.ix(j) = ix; 
                            suppf.iy(j) = iy;
                            suppf.iz(j) = iz;
                            j = j+1;
                        end
                    end                
                end
                % add the nodes in the interior y layers 
                if(nx ~= nxf)
                    for i2 = 1:n3
                        v1 = (nz-1)*(nxy1+n3*nxy2) + (nx-1)*(ny1+n3*ny2) + 1 + (ny1+(i2-1)*ny2) + (ny0-1);
                        for j2 = 0:(nyf-ny0)
                            suppf.nodes(j) = v1+j2; 
                            suppf.ix(j) = ix; 
                            suppf.iy(j) = iy;
                            suppf.iz(j) = iz;
                            j = j+1;
                        end
                    end
                end
            end
            % add the nodes in the interior xy layers
            if(nz ~= nzf)
                for i2 = 1:n3
                    v1 = (nz-1)*(nxy1+n3*nxy2) + 1 + (nxy1+(i2-1)*nxy2) + (nx0-1)*ny2 + (ny0-1);
                    for i3 = 0:(nxf-nx0)
                        v2 = v1+i3*ny2;
                        for j2 = 0:(nyf-ny0)
                            suppf.nodes(j) = v2+j2; 
                            suppf.ix(j) = ix; 
                            suppf.iy(j) = iy;
                            suppf.iz(j) = iz;
                            j = j+1;
                        end
                    end
                end
            end
        end
    end

    s = 3*suppf.nodes;
    ind1 = [s-2; s-1; s];
    ind1 = ind1(:);
    ind2 = [suppf.ix; suppf.iy; suppf.iz]; 
    ind2 = ind2(:);
    suppf.ind = ind1(ind2~=0); 
    str.freedofs = 1:3*nnodes; 
    str.freedofs(suppf.ind) = [];
    % ---------- %

    % Converting Loads % 
    j = 1;
    for i = 1:size(str.loads,1)
        % Initial and final coordinates of the load aplication
        x0 = str.loads(i,1);
        xf = str.loads(i,2);
        y0 = str.loads(i,3);
        yf = str.loads(i,4);
        z0 = str.loads(i,5);
        zf = str.loads(i,6);
        % Components of load
        Fx = str.loads(i,7); 
        Fy = str.loads(i,8); 
        Fz = str.loads(i,9);
        
        nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Initial nodes layer (on the linear grid) of load application in the x direction 
        nxf = 1 + floor(xf/str.el) + (mod(xf,str.el) > str.el/2); % Final nodes layer (on the linear grid) of load application in the x direction 
        ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Initial nodes layer (on the linear grid) of load application in the y direction
        nyf = 1 + floor(yf/str.eh) + (mod(yf,str.eh) > str.eh/2); % Final nodes layer (on the linear grid) of load application in the y direction
        nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Initial nodes layer (on the linear grid) of load application in the z direction
        nzf = 1 + floor(zf/str.ew) + (mod(zf,str.ew) > str.ew/2); % Final nodes layer (on the linear grid) of load application in the z direction        
        nx0s = 1 + floor(x0/elnew) + (mod(x0,elnew) > elnew/2); % Initial nodes layer of load application in the x direction 
        nxfs = 1 + floor(xf/elnew) + (mod(xf,elnew) > elnew/2); % Final nodes layer of load application in the x direction 
        ny0s = 1 + floor(y0/ehnew) + (mod(y0,ehnew) > ehnew/2); % Initial nodes layer of load application in the y direction
        nyfs = 1 + floor(yf/ehnew) + (mod(yf,ehnew) > ehnew/2); % Final nodes layer of load application in the y direction
        nz0s = 1 + floor(z0/ewnew) + (mod(z0,ewnew) > ewnew/2); % Initial nodes layer of load application in the z direction
        nzfs = 1 + floor(zf/ewnew) + (mod(zf,ewnew) > ewnew/2); % Final nodes layer of load application in the z direction

        node = (nz0-1)*(nxy1+n3*nxy2) + (nx0-1)*(ny1+n3*ny2) + (ny0-1)*(1+n3) + 1; % Initial node where the load is applied
        
        % Concentrated Load %
        if(x0 == xf && y0 == yf && z0 == zf)
            loads.nodes(j) = node;
            loads.Fx(j) = Fx;
            loads.Fy(j) = Fy;
            loads.Fz(j) = Fz;
            j = j+1; 
            continue;
        end
        
        % Load distributed only in the x direction %
        if(x0 ~= xf && y0 == yf && z0 == zf)
            if(nx0s == nxfs)
                loads.nodes(j) = node; 
                loads.Fx(j) = (xf-x0)*Fx;
                loads.Fy(j) = (xf-x0)*Fy;
                loads.Fz(j) = (xf-x0)*Fz;
                j = j+1;
            else
                loads.nodes(j) = node; 
                loads.Fx(j) = ((nx0s-1/2)*elnew - x0)*Fx;
                loads.Fy(j) = ((nx0s-1/2)*elnew - x0)*Fy;
                loads.Fz(j) = ((nx0s-1/2)*elnew - x0)*Fz;
                j = j+1;          
                for i2 = (nx0+1):(nxf-1)
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = elnew*Fx;
                        loads.Fy(j) = elnew*Fy;
                        loads.Fz(j) = elnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + ny1+n3*ny2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = elnew*Fx;
                    loads.Fy(j) = elnew*Fy;
                    loads.Fz(j) = elnew*Fz;
                    j = j+1;
                end
                node0 = node;
                for i3 = 1:n3
                    node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = elnew*Fx;
                    loads.Fy(j) = elnew*Fy;
                    loads.Fz(j) = elnew*Fz;
                    j = j+1; 
                end
                node = node0 + ny1+n3*ny2;
                loads.nodes(j) = node; 
                loads.Fx(j) = (xf - (nxfs-3/2)*elnew)*Fx;
                loads.Fy(j) = (xf - (nxfs-3/2)*elnew)*Fy;
                loads.Fz(j) = (xf - (nxfs-3/2)*elnew)*Fz;
                j = j+1;
            end  
            continue;
        end
        
        % Load distributed only in the y direction %
        if(x0 == xf && y0 ~= yf && z0 == zf)
            if(ny0s == nyfs)
                loads.nodes(j) = node; 
                loads.Fx(j) = (yf-y0)*Fx;
                loads.Fy(j) = (yf-y0)*Fy;
                loads.Fz(j) = (yf-y0)*Fz;
                j = j+1;
            else
                loads.nodes(j) = node; 
                loads.Fx(j) = ((ny0s-1/2)*ehnew - y0)*Fx;
                loads.Fy(j) = ((ny0s-1/2)*ehnew - y0)*Fy;
                loads.Fz(j) = ((ny0s-1/2)*ehnew - y0)*Fz;
                j = j+1;          
                for i2 = (ny0+1):(nyf-1)
                    for i3 = 1:n3
                        node = node + 1;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = ehnew*Fx;
                        loads.Fy(j) = ehnew*Fy;
                        loads.Fz(j) = ehnew*Fz;
                        j = j+1; 
                    end
                    node = node + 1;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = ehnew*Fx;
                    loads.Fy(j) = ehnew*Fy;
                    loads.Fz(j) = ehnew*Fz;
                    j = j+1;
                end
                for i3 = 1:n3
                    node = node + 1;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = ehnew*Fx;
                    loads.Fy(j) = ehnew*Fy;
                    loads.Fz(j) = ehnew*Fz;
                    j = j+1; 
                end
                node = node + 1;
                loads.nodes(j) = node; 
                loads.Fx(j) = (yf - (nyfs-3/2)*ehnew)*Fx;
                loads.Fy(j) = (yf - (nyfs-3/2)*ehnew)*Fy;
                loads.Fz(j) = (yf - (nyfs-3/2)*ehnew)*Fz;
                j = j+1;
            end 
            continue;
        end
        
        % Load distributed only in the z direction %
        if(x0 == xf && y0 == yf && z0 ~= zf)          
            if(nz0s == nzfs)
                loads.nodes(j) = node; 
                loads.Fx(j) = (zf-z0)*Fx;
                loads.Fy(j) = (zf-z0)*Fy;
                loads.Fz(j) = (zf-z0)*Fz;
                j = j+1;
            else
                loads.nodes(j) = node; 
                loads.Fx(j) = ((nz0s-1/2)*ewnew - z0)*Fx;
                loads.Fy(j) = ((nz0s-1/2)*ewnew - z0)*Fy;
                loads.Fz(j) = ((nz0s-1/2)*ewnew - z0)*Fz;
                j = j+1;          
                for i2 = (nz0+1):(nzf-1)
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = ewnew*Fx;
                        loads.Fy(j) = ewnew*Fy;
                        loads.Fz(j) = ewnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + nxy1+n3*nxy2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = ewnew*Fx;
                    loads.Fy(j) = ewnew*Fy;
                    loads.Fz(j) = ewnew*Fz;
                    j = j+1;
                end
                node0 = node;
                for i3 = 1:n3
                    node = node0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = ewnew*Fx;
                    loads.Fy(j) = ewnew*Fy;
                    loads.Fz(j) = ewnew*Fz;
                    j = j+1; 
                end
                node = node0 + nxy1+n3*nxy2;
                loads.nodes(j) = node; 
                loads.Fx(j) = (zf - (nzfs-3/2)*ewnew)*Fx;
                loads.Fy(j) = (zf - (nzfs-3/2)*ewnew)*Fy;
                loads.Fz(j) = (zf - (nzfs-3/2)*ewnew)*Fz;
                j = j+1;
            end     
            continue;
        end
        
        % Load distributed in an area parallel to the xy plane %
        if(x0 ~= xf && y0 ~= yf && z0 == zf)
            if(nx0s == nxfs)
                d = xf-x0;
                if(ny0s == nyfs)
                    loads.nodes(j) = node;
                    loads.Fx(j) = d*(yf-y0)*Fx; 
                    loads.Fy(j) = d*(yf-y0)*Fy;
                    loads.Fz(j) = d*(yf-y0)*Fz;
                    j = j+1; 
                else
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0s-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0s-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0s-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    for i2 = (ny0+1):(nyf-1)
                        for i3 = 1:n3
                            node = node + 1;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ehnew*Fx;
                            loads.Fy(j) = d*ehnew*Fy;
                            loads.Fz(j) = d*ehnew*Fz;
                            j = j+1; 
                        end
                        node = node + 1;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1;
                    end
                    for i3 = 1:n3
                        node = node + 1;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                    end
                    node = node + 1;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyfs-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyfs-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyfs-3/2)*ehnew)*Fz;
                    j = j+1;
                end 
            else
                if(ny0s == nyfs)
                    d = yf-y0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0s-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0s-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0s-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    for i2 = (nx0+1):(nxf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + ny1+n3*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + ny1+n3*ny2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxfs-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxfs-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxfs-3/2)*elnew)*Fz;
                    j = j+1;
                else   
                    n0 = node; 
                    d = (ny0s-1/2)*ehnew - y0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0s-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0s-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0s-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    for i2 = (nx0+1):(nxf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + ny1+n3*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + ny1+n3*ny2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxf-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxf-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxf-3/2)*elnew)*Fz;
                    j = j+1;

                    d = ehnew;
                    for j2 = (ny0+1):(nyf-1)
                        for i3 = 1:n3
                            node = n0 + i3;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*((nx0-1/2)*str.el - x0)*Fx;
                            loads.Fy(j) = d*((nx0-1/2)*str.el - x0)*Fy;
                            loads.Fz(j) = d*((nx0-1/2)*str.el - x0)*Fz;
                            j = j+1;
                            node = node + ny2; 
                            for i2 = (nx0+1):(nxf-1)
                                loads.nodes(j) = node; 
                                loads.Fx(j) = d*str.el*Fx;
                                loads.Fy(j) = d*str.el*Fy;
                                loads.Fz(j) = d*str.el*Fz;
                                j = j+1; 
                                node = node + ny2;
                            end
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*(xf - (nxf-3/2)*str.el)*Fx;
                            loads.Fy(j) = d*(xf - (nxf-3/2)*str.el)*Fy;
                            loads.Fz(j) = d*(xf - (nxf-3/2)*str.el)*Fz;
                            j = j+1;
                        end
                        n0 = n0 + elemDeg;
                        node = n0;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*((nx0s-1/2)*elnew - x0)*Fx;
                        loads.Fy(j) = d*((nx0s-1/2)*elnew - x0)*Fy;
                        loads.Fz(j) = d*((nx0s-1/2)*elnew - x0)*Fz;
                        j = j+1; 
                        for i2 = (nx0+1):(nxf-1)
                            node0 = node;
                            for i3 = 1:n3
                                node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                                loads.nodes(j) = node; 
                                loads.Fx(j) = d*elnew*Fx;
                                loads.Fy(j) = d*elnew*Fy;
                                loads.Fz(j) = d*elnew*Fz;
                                j = j+1; 
                            end
                            node = node0 + ny1+n3*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + ny1+n3*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*(xf - (nxfs-3/2)*elnew)*Fx;
                        loads.Fy(j) = d*(xf - (nxfs-3/2)*elnew)*Fy;
                        loads.Fz(j) = d*(xf - (nxfs-3/2)*elnew)*Fz;
                        j = j+1;
                    end
                    
                    for i3 = 1:n3
                        node = n0 + i3;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*((nx0-1/2)*str.el - x0)*Fx;
                        loads.Fy(j) = d*((nx0-1/2)*str.el - x0)*Fy;
                        loads.Fz(j) = d*((nx0-1/2)*str.el - x0)*Fz;
                        j = j+1;
                        node = node + ny2; 
                        for i2 = (nx0+1):(nxf-1)
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*str.el*Fx;
                            loads.Fy(j) = d*str.el*Fy;
                            loads.Fz(j) = d*str.el*Fz;
                            j = j+1; 
                            node = node + ny2;
                        end
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*(xf - (nxf-3/2)*str.el)*Fx;
                        loads.Fy(j) = d*(xf - (nxf-3/2)*str.el)*Fy;
                        loads.Fz(j) = d*(xf - (nxf-3/2)*str.el)*Fz;
                        j = j+1;
                    end
                    node = n0 + elemDeg;
                    d = yf - (nyfs-3/2)*ehnew;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0s-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0s-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0s-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    for i2 = (nx0+1):(nxf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + ny1+n3*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + ny1+n3*ny2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxfs-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxfs-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxfs-3/2)*elnew)*Fz;
                    j = j+1;
                end 
            end
            continue;
        end
        
        % Load distributed in an area parallel to the xz plane %
        if(x0 ~= xf && y0 == yf && z0 ~= zf)
            if(nx0s == nxfs)
                d = xf-x0;
                if(nz0s == nzfs)
                    loads.nodes(j) = node;
                    loads.Fx(j) = d*(zf-z0)*Fx; 
                    loads.Fy(j) = d*(zf-z0)*Fy;
                    loads.Fz(j) = d*(zf-z0)*Fz;
                    j = j+1; 
                else
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nz0s-1/2)*ewnew - z0)*Fx;
                    loads.Fy(j) = d*((nz0s-1/2)*ewnew - z0)*Fy;
                    loads.Fz(j) = d*((nz0s-1/2)*ewnew - z0)*Fz;
                    j = j+1;
                    for i2 = (nz0+1):(nzf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ewnew*Fx;
                            loads.Fy(j) = d*ewnew*Fy;
                            loads.Fz(j) = d*ewnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + nxy1+n3*nxy2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ewnew*Fx;
                        loads.Fy(j) = d*ewnew*Fy;
                        loads.Fz(j) = d*ewnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ewnew*Fx;
                        loads.Fy(j) = d*ewnew*Fy;
                        loads.Fz(j) = d*ewnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + nxy1+n3*nxy2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(zf - (nzfs-3/2)*ewnew)*Fx;
                    loads.Fy(j) = d*(zf - (nzfs-3/2)*ewnew)*Fy;
                    loads.Fz(j) = d*(zf - (nzfs-3/2)*ewnew)*Fz;
                    j = j+1;
                end 
            else
                if(nz0s == nzfs)
                    d = zf-z0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0s-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0s-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0s-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    for i2 = (nx0+1):(nxf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + ny1+n3*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + ny1+n3*ny2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxfs-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxfs-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxfs-3/2)*elnew)*Fz;
                    j = j+1;
                else   
                    n0 = node; 
                    d = (nz0s-1/2)*ewnew - z0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0s-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0s-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0s-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    for i2 = (nx0+1):(nxf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + ny1+n3*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + ny1+n3*ny2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxf-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxf-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxf-3/2)*elnew)*Fz;
                    j = j+1;

                    d = ewnew;
                    for j2 = (nz0+1):(nzf-1)
                        for i3 = 1:n3
                            node = n0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*((nx0-1/2)*str.el - x0)*Fx;
                            loads.Fy(j) = d*((nx0-1/2)*str.el - x0)*Fy;
                            loads.Fz(j) = d*((nx0-1/2)*str.el - x0)*Fz;
                            j = j+1;
                            node = node + ny2; 
                            for i2 = (nx0+1):(nxf-1)
                                loads.nodes(j) = node; 
                                loads.Fx(j) = d*str.el*Fx;
                                loads.Fy(j) = d*str.el*Fy;
                                loads.Fz(j) = d*str.el*Fz;
                                j = j+1; 
                                node = node + ny2;
                            end
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*(xf - (nxf-3/2)*str.el)*Fx;
                            loads.Fy(j) = d*(xf - (nxf-3/2)*str.el)*Fy;
                            loads.Fz(j) = d*(xf - (nxf-3/2)*str.el)*Fz;
                            j = j+1;
                        end
                        n0 = n0 + nxy1+n3*nxy2;
                        node = n0;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*((nx0s-1/2)*elnew - x0)*Fx;
                        loads.Fy(j) = d*((nx0s-1/2)*elnew - x0)*Fy;
                        loads.Fz(j) = d*((nx0s-1/2)*elnew - x0)*Fz;
                        j = j+1; 
                        for i2 = (nx0+1):(nxf-1)
                            node0 = node;
                            for i3 = 1:n3
                                node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                                loads.nodes(j) = node; 
                                loads.Fx(j) = d*elnew*Fx;
                                loads.Fy(j) = d*elnew*Fy;
                                loads.Fz(j) = d*elnew*Fz;
                                j = j+1; 
                            end
                            node = node0 + ny1+n3*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + ny1+n3*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*(xf - (nxfs-3/2)*elnew)*Fx;
                        loads.Fy(j) = d*(xf - (nxfs-3/2)*elnew)*Fy;
                        loads.Fz(j) = d*(xf - (nxfs-3/2)*elnew)*Fz;
                        j = j+1;
                    end
                    
                    for i3 = 1:n3
                        node = n0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*((nx0-1/2)*str.el - x0)*Fx;
                        loads.Fy(j) = d*((nx0-1/2)*str.el - x0)*Fy;
                        loads.Fz(j) = d*((nx0-1/2)*str.el - x0)*Fz;
                        j = j+1;
                        node = node + ny2; 
                        for i2 = (nx0+1):(nxf-1)
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*str.el*Fx;
                            loads.Fy(j) = d*str.el*Fy;
                            loads.Fz(j) = d*str.el*Fz;
                            j = j+1; 
                            node = node + ny2;
                        end
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*(xf - (nxf-3/2)*str.el)*Fx;
                        loads.Fy(j) = d*(xf - (nxf-3/2)*str.el)*Fy;
                        loads.Fz(j) = d*(xf - (nxf-3/2)*str.el)*Fz;
                        j = j+1;
                    end
                    node = n0 + nxy1+n3*nxy2;
                    d = zf - (nzfs-3/2)*ewnew;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0s-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0s-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0s-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    for i2 = (nx0+1):(nxf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*elnew*Fx;
                            loads.Fy(j) = d*elnew*Fy;
                            loads.Fz(j) = d*elnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + ny1+n3*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + ((str.nely-ny0+1)*elemDeg+ny0) + (i3-1)*ny2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + ny1+n3*ny2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxfs-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxfs-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxfs-3/2)*elnew)*Fz;
                    j = j+1;
                end 
            end
            continue;
        end
        
        % Load distributed in an area parallel to the yz plane %
        if(x0 == xf && y0 ~= yf && z0 ~= zf)
            if(ny0s == nyfs)
                d = yf-y0;
                if(nz0s == nzfs)
                    loads.nodes(j) = node;
                    loads.Fx(j) = d*(zf-z0)*Fx; 
                    loads.Fy(j) = d*(zf-z0)*Fy;
                    loads.Fz(j) = d*(zf-z0)*Fz;
                    j = j+1; 
                else
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nz0s-1/2)*ewnew - z0)*Fx;
                    loads.Fy(j) = d*((nz0s-1/2)*ewnew - z0)*Fy;
                    loads.Fz(j) = d*((nz0s-1/2)*ewnew - z0)*Fz;
                    j = j+1;
                    for i2 = (nz0+1):(nzf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ewnew*Fx;
                            loads.Fy(j) = d*ewnew*Fy;
                            loads.Fz(j) = d*ewnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + nxy1+n3*nxy2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ewnew*Fx;
                        loads.Fy(j) = d*ewnew*Fy;
                        loads.Fz(j) = d*ewnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ewnew*Fx;
                        loads.Fy(j) = d*ewnew*Fy;
                        loads.Fz(j) = d*ewnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + nxy1+n3*nxy2;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(zf - (nzfs-3/2)*ewnew)*Fx;
                    loads.Fy(j) = d*(zf - (nzfs-3/2)*ewnew)*Fy;
                    loads.Fz(j) = d*(zf - (nzfs-3/2)*ewnew)*Fz;
                    j = j+1;
                end 
            else
                if(nz0s == nzfs)
                    d = zf-z0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0s-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0s-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0s-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    for i2 = (ny0+1):(nyf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + i3;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ehnew*Fx;
                            loads.Fy(j) = d*ehnew*Fy;
                            loads.Fz(j) = d*ehnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + elemDeg;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + i3;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + elemDeg;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyfs-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyfs-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyfs-3/2)*ehnew)*Fz;
                    j = j+1;
                else   
                    n0 = node; 
                    d = (nz0s-1/2)*ewnew - z0;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0s-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0s-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0s-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    for i2 = (ny0+1):(nyf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + i3;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ehnew*Fx;
                            loads.Fy(j) = d*ehnew*Fy;
                            loads.Fz(j) = d*ehnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + elemDeg;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + i3;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + elemDeg;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyf-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyf-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyf-3/2)*ehnew)*Fz;
                    j = j+1;

                    d = ewnew;
                    for j2 = (nz0+1):(nzf-1)
                        for i3 = 1:n3
                            node = n0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*((ny0-1/2)*str.eh - y0)*Fx;
                            loads.Fy(j) = d*((ny0-1/2)*str.eh - y0)*Fy;
                            loads.Fz(j) = d*((ny0-1/2)*str.eh - y0)*Fz;
                            j = j+1;
                            node = node + 1; 
                            for i2 = (ny0+1):(nyf-1)
                                loads.nodes(j) = node; 
                                loads.Fx(j) = d*str.eh*Fx;
                                loads.Fy(j) = d*str.eh*Fy;
                                loads.Fz(j) = d*str.eh*Fz;
                                j = j+1; 
                                node = node + 1;
                            end
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                            loads.Fy(j) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                            loads.Fz(j) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                            j = j+1;
                        end
                        n0 = n0 + nxy1+n3*nxy2;
                        node = n0;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*((ny0s-1/2)*ehnew - y0)*Fx;
                        loads.Fy(j) = d*((ny0s-1/2)*ehnew - y0)*Fy;
                        loads.Fz(j) = d*((ny0s-1/2)*ehnew - y0)*Fz;
                        j = j+1; 
                        for i2 = (ny0+1):(nyf-1)
                            node0 = node;
                            for i3 = 1:n3
                                node = node0 + i3;
                                loads.nodes(j) = node; 
                                loads.Fx(j) = d*ehnew*Fx;
                                loads.Fy(j) = d*ehnew*Fy;
                                loads.Fz(j) = d*ehnew*Fz;
                                j = j+1; 
                            end
                            node = node0 + elemDeg;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ehnew*Fx;
                            loads.Fy(j) = d*ehnew*Fy;
                            loads.Fz(j) = d*ehnew*Fz;
                            j = j+1; 
                        end
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + i3;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ehnew*Fx;
                            loads.Fy(j) = d*ehnew*Fy;
                            loads.Fz(j) = d*ehnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + elemDeg;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*(yf - (nyfs-3/2)*ehnew)*Fx;
                        loads.Fy(j) = d*(yf - (nyfs-3/2)*ehnew)*Fy;
                        loads.Fz(j) = d*(yf - (nyfs-3/2)*ehnew)*Fz;
                        j = j+1;
                    end
                    
                    for i3 = 1:n3
                        node = n0 + (nxy1-((nx0-1)*(ny1+n3*ny2)) + (nx0-1)*ny2 + (ny0-1)*(1-elemDeg)) + (i3-1)*nxy2;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*((ny0-1/2)*str.eh - y0)*Fx;
                        loads.Fy(j) = d*((ny0-1/2)*str.eh - y0)*Fy;
                        loads.Fz(j) = d*((ny0-1/2)*str.eh - y0)*Fz;
                        j = j+1;
                        node = node + 1; 
                        for i2 = (ny0+1):(nyf-1)
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*str.eh*Fx;
                            loads.Fy(j) = d*str.eh*Fy;
                            loads.Fz(j) = d*str.eh*Fz;
                            j = j+1; 
                            node = node + 1;
                        end
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                        loads.Fy(j) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                        loads.Fz(j) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                        j = j+1;
                    end
                    node = n0 + nxy1+n3*nxy2;
                    d = zf - (nzfs-3/2)*ewnew;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((ny0s-1/2)*ehnew - y0)*Fx;
                    loads.Fy(j) = d*((ny0s-1/2)*ehnew - y0)*Fy;
                    loads.Fz(j) = d*((ny0s-1/2)*ehnew - y0)*Fz;
                    j = j+1;
                    for i2 = (ny0+1):(nyf-1)
                        node0 = node;
                        for i3 = 1:n3
                            node = node0 + i3;
                            loads.nodes(j) = node; 
                            loads.Fx(j) = d*ehnew*Fx;
                            loads.Fy(j) = d*ehnew*Fy;
                            loads.Fz(j) = d*ehnew*Fz;
                            j = j+1; 
                        end
                        node = node0 + elemDeg;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1;
                    end
                    node0 = node;
                    for i3 = 1:n3
                        node = node0 + i3;
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*ehnew*Fx;
                        loads.Fy(j) = d*ehnew*Fy;
                        loads.Fz(j) = d*ehnew*Fz;
                        j = j+1; 
                    end
                    node = node0 + elemDeg;
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(yf - (nyfs-3/2)*ehnew)*Fx;
                    loads.Fy(j) = d*(yf - (nyfs-3/2)*ehnew)*Fy;
                    loads.Fz(j) = d*(yf - (nyfs-3/2)*ehnew)*Fz;
                    j = j+1;
                end 
            end
            continue;
        end
    end

    str.f = zeros(3*nnodes, 1); % Vector of nodal loads
    s = 3*loads.nodes; 
    for i = 1:length(s)
        str.f(s(i)-2) = str.f(s(i)-2) + loads.Fx(i); 
        str.f(s(i)-1) = str.f(s(i)-1) + loads.Fy(i);
        str.f(s(i)) = str.f(s(i)) + loads.Fz(i);
    end 
    str.f(suppf.ind) = []; % Boundary conditions
    % ---------- %
    str.f = str.f/norm(str.f);

    str.nnodes = nnodes; % Adapts the total number of nodes on the grid
end

end