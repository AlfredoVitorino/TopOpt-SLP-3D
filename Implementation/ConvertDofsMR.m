function [strT] = ConvertDofsMR(str,mr)
% ConvertDofsMR adapts the number of nodes, freedofs and the vector of nodal loads to the density mesh used in multiresolution. %
% INPUT: str - structure with the problem data.
%        elem - structure with element characteristics (deg - polynomial degree, type - element type (1 = Lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: strT - scruture with the problem data converted to the density mesh.
% ---------- %

strT.l = str.l*(mr.n);
strT.h = str.h*(mr.n);
strT.w = str.w*(mr.n);
strT.nelx = str.nelx*(mr.n); 
strT.nely = str.nely*(mr.n); 
strT.nelz = str.nelz*(mr.n);
strT.el = strT.l/strT.nelx; 
strT.eh = strT.h/strT.nely; 
strT.ew = strT.w/strT.nelz; 
strT.nx = strT.nelx+1; 
strT.ny = strT.nely+1;
strT.nz = strT.nelz+1;
strT.nelem = strT.nelx*strT.nely*strT.nelz;
strT.nnodes = strT.nx*strT.ny*strT.nz; 
strT.E = str.E; 
strT.nu = str.nu;
elnew = str.el/(mr.n); % distance between two nodes in the x direction 
ehnew = str.eh/(mr.n); % distance between two nodes in the y direction 
ewnew = str.ew/(mr.n); % distance between two nodes in the z direction 

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

    node = (nz0-1)*strT.nx*strT.ny + (nx0-1)*strT.ny + ny0; % Initial node of the support
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
            node = node - 1 + (strT.ny - nyf) + ny0;
        end
        node = node + strT.ny*(strT.nx - nxf + nx0 - 1);
    end
end
s = 3*supp.nodes;
ind1 = [s-2; s-1; s];
ind1 = ind1(:);
ind2 = [supp.ix; supp.iy; supp.iz]; 
ind2 = ind2(:);
supp.ind = ind1(ind2~=0); % Global indexes of the nodes with supports preventing the movement
strT.freedofs = 1:3*strT.nnodes; 
strT.freedofs(supp.ind) = []; % Indexes of nonzero degrees of freedom. 
% ---------- %

% Converting loads %
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

    nx0 = 1 + floor(x0/elnew) + (mod(x0,elnew) > elnew/2); % Initial nodes layer of load application in the x direction 
    nxf = 1 + floor(xf/elnew) + (mod(xf,elnew) > elnew/2); % Final nodes layer of load application in the x direction 
    ny0 = 1 + floor(y0/ehnew) + (mod(y0,ehnew) > ehnew/2); % Initial nodes layer of load application in the y direction
    nyf = 1 + floor(yf/ehnew) + (mod(yf,ehnew) > ehnew/2); % Final nodes layer of load application in the y direction
    nz0 = 1 + floor(z0/ewnew) + (mod(z0,ewnew) > ewnew/2); % Initial nodes layer of load application in the z direction
    nzf = 1 + floor(zf/ewnew) + (mod(zf,ewnew) > ewnew/2); % Final nodes layer of load application in the z direction

    node = (nz0-1)*strT.nx*strT.ny + (nx0-1)*strT.ny + ny0; % Initial node where the load is applied

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
            node = node + strT.ny; 
            for i2 = (nx0+1):(nxf-1)
                loads.nodes(j) = node; 
                loads.Fx(j) = elnew*Fx;
                loads.Fy(j) = elnew*Fy;
                loads.Fz(j) = elnew*Fz;
                j = j+1; 
                node = node + strT.ny;
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
            node = node + strT.nx*strT.ny; 
            for i2 = (nz0+1):(nzf-1)
                loads.nodes(j) = node; 
                loads.Fx(j) = ewnew*Fx;
                loads.Fy(j) = ewnew*Fy;
                loads.Fz(j) = ewnew*Fz;
                j = j+1; 
                node = node + strT.nx*strT.ny;
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
                node = node + strT.ny; 
                for i2 = (nx0+1):(nxf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*elnew*Fx;
                    loads.Fy(j) = d*elnew*Fy;
                    loads.Fz(j) = d*elnew*Fz;
                    j = j+1; 
                    node = node + strT.ny;
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
                node = node + (strT.ny - nyf) + ny0;

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
                    node = node + (strT.ny - nyf) + ny0;
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
                node = node + strT.nx*strT.ny; 
                for i2 = (nz0+1):(nzf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*ewnew*Fx;
                    loads.Fy(j) = d*ewnew*Fy;
                    loads.Fz(j) = d*ewnew*Fz;
                    j = j+1; 
                    node = node + strT.nx*strT.ny;
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
                node = node + strT.ny; 
                for i2 = (nx0+1):(nxf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*elnew*Fx;
                    loads.Fy(j) = d*elnew*Fy;
                    loads.Fz(j) = d*elnew*Fz;
                    j = j+1; 
                    node = node + strT.ny;
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
                node = node + strT.ny; 
                for i2 = (nx0+1):(nxf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*elnew*Fx;
                    loads.Fy(j) = d*elnew*Fy;
                    loads.Fz(j) = d*elnew*Fz;
                    j = j+1; 
                    node = node + strT.ny;
                end
                loads.nodes(j) = node; 
                loads.Fx(j) = d*(xf - (nxf-3/2)*elnew)*Fx;
                loads.Fy(j) = d*(xf - (nxf-3/2)*elnew)*Fy;
                loads.Fz(j) = d*(xf - (nxf-3/2)*elnew)*Fz;
                j = j+1;
                node = node + strT.nx*strT.ny - (nxf-nx0)*strT.ny;

                d = ewnew;
                for j2 = (nz0+1):(nzf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*((nx0-1/2)*elnew - x0)*Fx;
                    loads.Fy(j) = d*((nx0-1/2)*elnew - x0)*Fy;
                    loads.Fz(j) = d*((nx0-1/2)*elnew - x0)*Fz;
                    j = j+1;
                    node = node + strT.ny; 
                    for i2 = (nx0+1):(nxf-1)
                        loads.nodes(j) = node; 
                        loads.Fx(j) = d*elnew*Fx;
                        loads.Fy(j) = d*elnew*Fy;
                        loads.Fz(j) = d*elnew*Fz;
                        j = j+1; 
                        node = node + strT.ny;
                    end
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*(xf - (nxf-3/2)*elnew)*Fx;
                    loads.Fy(j) = d*(xf - (nxf-3/2)*elnew)*Fy;
                    loads.Fz(j) = d*(xf - (nxf-3/2)*elnew)*Fz;
                    j = j+1;
                    node = node + strT.nx*strT.ny - (nxf-nx0)*strT.ny;
                end

                d = zf - (nzf-3/2)*ewnew;
                loads.nodes(j) = node; 
                loads.Fx(j) = d*((nx0-1/2)*elnew - x0)*Fx;
                loads.Fy(j) = d*((nx0-1/2)*elnew - x0)*Fy;
                loads.Fz(j) = d*((nx0-1/2)*elnew - x0)*Fz;
                j = j+1;
                node = node + strT.ny; 
                for i2 = (nx0+1):(nxf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*elnew*Fx;
                    loads.Fy(j) = d*elnew*Fy;
                    loads.Fz(j) = d*elnew*Fz;
                    j = j+1; 
                    node = node + strT.ny;
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
                node = node + strT.nx*strT.ny; 
                for i2 = (nz0+1):(nzf-1)
                    loads.nodes(j) = node; 
                    loads.Fx(j) = d*ewnew*Fx;
                    loads.Fy(j) = d*ewnew*Fy;
                    loads.Fz(j) = d*ewnew*Fz;
                    j = j+1; 
                    node = node + strT.nx*strT.ny;
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
                node = node + strT.nx*strT.ny - (nyf-ny0);

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
                    node = node + strT.nx*strT.ny - (nyf-ny0);
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

strT.f = zeros(3*strT.nnodes, 1); % Vector of nodal loads
s = 3*loads.nodes; 
for i = 1:length(s)
    strT.f(s(i)-2) = strT.f(s(i)-2) + loads.Fx(i); 
    strT.f(s(i)-1) = strT.f(s(i)-1) + loads.Fy(i);
    strT.f(s(i)) = strT.f(s(i)) + loads.Fz(i);
end 
strT.f(supp.ind) = []; % Boundary conditions
% ---------- %

end