function [str] = EnterData
% EnterData reads and constructs the problem data. %
% OUTPUT: str - structure with the problem data.
% ---------- % 

% Domain dimensions %
while(1)
    str.l = input('Length: '); 
    if(isempty(str.l) || str.l <= 0)
        fprintf('Error: the length must be positive. \n'); 
    else
        break;
    end
end
while(1)
    str.h = input('Height: ');
    if(isempty(str.h) || str.h <= 0)
        fprintf('Error: the height must be positive. \n'); 
    else
        break;
    end       
end
while(1)
    str.w = input('Width: ');
    if(isempty(str.w) || str.w <= 0)
        fprintf('Error: the width must be positive. \n');
    else
        break;
    end
end
% ---------- %

% Number of finite elements % 
while(1)
    str.nelx = input('Number of elements in the length(x) direction: ');
    if(isempty(str.nelx) || str.nelx <= 0)
        fprintf('Error: this value must be positive. \n');
    else
        break;
    end
end
while(1)
    str.nely = input('Number of elements in the height(y) direction: ');
    if(isempty(str.nely) || str.nely <= 0)
        fprintf('Error: this value must be positive. \n');
    else
        break;
    end
end
while(1)
    str.nelz = input('Number of elements in the width(z) direction: ');
    if(isempty(str.nelz) || str.nelz <= 0)
        fprintf('Error: this value must be positive. \n');
    else
        break;
    end
end
% ---------- %

% Material parameters %
while(1)
    str.E = input('Young`s modulus: ');
    if(isempty(str.E) || str.E <= 0)
        fprintf('Error: this value must be positive. \n');
    else
        break;
    end
end
while(1)
    str.nu = input('Poisson`s ratio: ');
    if(isempty(str.nu) || str.nu <= 0 || str.nu > 1)
        fprintf('Error: this value must be between 0 and 1. \n');
    else
        break;
    end
end
% ---------- %

str.el = str.l/str.nelx; % Element length
str.eh = str.h/str.nely; % Element height
str.ew = str.w/str.nelz; % Element width
str.nx = str.nelx + 1; % Number of nodes in the length direction
str.ny = str.nely + 1; % Number of nodes in the height direction
str.nz = str.nelz + 1; % Number of nodes in the width direction
str.nelem = str.nelx*str.nely*str.nelz; % Total number of elements
str.nnodes = str.nx*str.ny*str.nz; % Total number of nodes

% Supports (boundary conditions) % 
p = 1;
str.supp = [];
while(1)
    fprintf('\n');
    op = input('Continue adding supports? (1 = Yes, 0 = No) ');
    if(op == 0)
        break;
    elseif(isempty(op) || op ~= 1)
        continue; 
    end
    
    x0 = input('Initial x-coordinate of the support: '); 
    xf = input('Final x-coordinate of the support: ');
    if((x0 < 0 && xf < 0) || (x0 > str.l && xf > str.l))
        fprintf('Error: support is outside the domain.\n');
        continue;
    end
    y0 = input('Initial y-coordinate of the support: '); 
    yf = input('Final y-coordinate of the support: '); 
    if((y0 < 0 && yf < 0) || (y0 > str.h && yf > str.h))
        fprintf('Error: support is outside the domain.\n');
        continue;
    end
    z0 = input('Initial z-coordinate of the support: '); 
    zf = input('Final z-coordinate of the support: '); 
    if((z0 < 0 && zf < 0) || (z0 > str.w && zf > str.w))
        fprintf('Error: support is outside the domain.\n');
        continue;
    end
    
    ix = input('Does the support prevent the movement in the x direction? (1 = Yes, 0 = No) ');
    iy = input('Does the support prevent the movement in the y direction? (1 = Yes, 0 = No) ');
    iz = input('Does the support prevent the movement in the z direction? (1 = Yes, 0 = No) ');
    
    nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Initial nodes layer of the support in the x direction 
    nxf = 1 + floor(xf/str.el) + (mod(xf,str.el) > str.el/2); % Final nodes layer of the support in the x direction
    ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Initial nodes layer of the support in the y direction
    nyf = 1 + floor(yf/str.eh) + (mod(yf,str.eh) > str.eh/2); % Final nodes layer of the support in the y direction
    nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Initial nodes layer of the support in the z direction
    nzf = 1 + floor(zf/str.ew) + (mod(zf,str.ew) > str.ew/2); % Final nodes layer of the support in the z direction
    
    node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Initial node of the support
    for i = nz0:nzf
        for j = nx0:nxf
            for k = ny0:nyf
                supp.nodes(p) = node;
                supp.ix(p) = ix;
                supp.iy(p) = iy;
                supp.iz(p) = iz;
                node = node + 1;
                p = p + 1;
            end
            node = node - 1 + (str.ny - nyf) + ny0;
        end
        node = node + str.ny*(str.nx - nxf + nx0 - 1);
    end
    
    str.supp = [str.supp; [x0,xf,y0,yf,z0,zf,ix,iy,iz]]; % Save informations about each support 
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

% Loads %
p = 1;
str.loads = [];
while(1)
    fprintf('\n');
    op = input('Continue adding loads? (1 = Yes, 0 = No) ');
    if(op == 0)
        break;
    elseif(isempty(op) || op ~= 1)
        continue; 
    end
    
    Fx = input('x-component of load: '); 
    Fy = input('y-component of load: '); 
    Fz = input('z-component of load: ');   
    opx = input('Is the load distributed in the x direction? (1 = Yes, 0 = No) ');
    opy = input('Is the load distributed in the y direction? (1 = Yes, 0 = No) ');
    opz = input('Is the load distributed in the z direction? (1 = Yes, 0 = No) ');
    
    % Concentrated load %
    if(opx == 0 && opy == 0 && opz == 0) 
        x0 = input('x-coordinate of load application: ');
        y0 = input('y-coordinate of load application: ');
        z0 = input('z-coordinate of load application: ');
        if(x0 < 0 || y0 < 0 || z0 < 0 || x0 > str.l || y0 > str.h || z0 > str.w) 
            fprintf('Error: the load is applied outside the domain.');
            continue;
        end
        
        nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Nodes layer of load application in the x direction
        ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Nodes layer of load application in the y direction
        nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Nodes layer of load application in the z direction
        node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Node of load application
        
        loads.nodes(p) = node;
        loads.Fx(p) = Fx;
        loads.Fy(p) = Fy;
        loads.Fz(p) = Fz;
        p = p+1; 
        
        str.loads = [str.loads; [x0,x0,y0,y0,z0,z0,Fx,Fy,Fz]]; % Save informations about the load
        continue; 
    end 
    
    % Load distributed only in the x direction %
    if(opx == 1 && opy == 0 && opz == 0) 
        x0 = input('Initial x-coordinate of load application: ');
        xf = input('Final x-coordinate of load application: ');
        y0 = input('y-coordinate of load application: ');
        z0 = input('z-coordinate of load application: ');
        if(x0 < 0 || xf < x0 || y0 < 0 || z0 < 0 || x0 > str.l || xf > str.l || y0 > str.h || z0 > str.w) 
            fprintf('Error: the load is applied outside the domain.');
            continue;
        end
        
        nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Initial nodes layer of load application in the x direction
        nxf = 1 + floor(xf/str.el) + (mod(xf,str.el) > str.el/2); % Final nodes layer of load application in the x direction
        ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Nodes layer of load application in the y direction
        nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Nodes layer of load application in the z direction
        node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Initial node of load application
        
        if(nx0 == nxf)
            loads.nodes(p) = node; 
            loads.Fx(p) = (xf-x0)*Fx;
            loads.Fy(p) = (xf-x0)*Fy;
            loads.Fz(p) = (xf-x0)*Fz;
            p = p+1;
        else
            loads.nodes(p) = node; 
            loads.Fx(p) = ((nx0-1/2)*str.el - x0)*Fx;
            loads.Fy(p) = ((nx0-1/2)*str.el - x0)*Fy;
            loads.Fz(p) = ((nx0-1/2)*str.el - x0)*Fz;
            p = p+1;
            node = node + str.ny; 
            for i = (nx0+1):(nxf-1)
                loads.nodes(p) = node; 
                loads.Fx(p) = str.el*Fx;
                loads.Fy(p) = str.el*Fy;
                loads.Fz(p) = str.el*Fz;
                p = p+1; 
                node = node + str.ny;
            end
            loads.nodes(p) = node; 
            loads.Fx(p) = (xf - (nxf-3/2)*str.el)*Fx;
            loads.Fy(p) = (xf - (nxf-3/2)*str.el)*Fy;
            loads.Fz(p) = (xf - (nxf-3/2)*str.el)*Fz;
            p = p+1;
        end
        
        str.loads = [str.loads; [x0,xf,y0,y0,z0,z0,Fx,Fy,Fz]]; % Save informations about the load
        continue;
    end 
    
    % Load distributed only in the y direction %
    if(opx == 0 && opy == 1 && opz == 0)
        x0 = input('x-coordinate of load application: ');
        y0 = input('Initial y-coordinate of load application: ');
        yf = input('Final y-coordinate of load application: ');
        z0 = input('z-coordinate of load application: ');
        if(y0 < 0 || yf < y0 || x0 < 0 || z0 < 0 || y0 > str.h || yf > str.h || x0 > str.l || z0 > str.w) 
            fprintf('Error: the load is applied outside the domain.');
            continue;
        end
        
        nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Nodes layer of load application in the x direction
        ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Initial nodes layer of load application in the y direction
        nyf = 1 + floor(yf/str.eh) + (mod(yf,str.eh) > str.eh/2); % Final nodes layer of load application in the y direction    
        nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Nodes layer of load application in the z direction
        node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Initial node of load application
        
        if(ny0 == nyf)
            loads.nodes(p) = node; 
            loads.Fx(p) = (yf-y0)*Fx;
            loads.Fy(p) = (yf-y0)*Fy;
            loads.Fz(p) = (yf-y0)*Fz;
            p = p+1;
        else
            loads.nodes(p) = node; 
            loads.Fx(p) = ((ny0-1/2)*str.eh - y0)*Fx;
            loads.Fy(p) = ((ny0-1/2)*str.eh - y0)*Fy;
            loads.Fz(p) = ((ny0-1/2)*str.eh - y0)*Fz;
            p = p+1;
            node = node + 1; 
            for i = (ny0+1):(nyf-1)
                loads.nodes(p) = node; 
                loads.Fx(p) = str.eh*Fx;
                loads.Fy(p) = str.eh*Fy;
                loads.Fz(p) = str.eh*Fz;
                p = p+1; 
                node = node + 1;
            end
            loads.nodes(p) = node; 
            loads.Fx(p) = (yf - (nyf-3/2)*str.eh)*Fx;
            loads.Fy(p) = (yf - (nyf-3/2)*str.eh)*Fy;
            loads.Fz(p) = (yf - (nyf-3/2)*str.eh)*Fz;
            p = p+1;
        end
        
        str.loads = [str.loads; [x0,x0,y0,yf,z0,z0,Fx,Fy,Fz]]; % Save informations about the load
        continue;
    end
    
    % Load distributed only in the z direction % 
    if(opx == 0 && opy == 0 && opz == 1) 
        x0 = input('x-coordinate of load application: ');
        y0 = input('y-coordinate of load application: ');
        z0 = input('Initial z-coordinate of load application: ');
        zf = input('Final z-coordinate of load application: ');
        if(z0 < 0 || zf < z0 || x0 < 0 || y0 < 0 || z0 > str.w || zf > str.w || x0 > str.l || y0 > str.h) 
            fprintf('Error: the load is applied outside the domain.');
            continue;
        end
        
        nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Nodes layer of load application in the x direction
        ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Nodes layer of load application in the y direction
        nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Initial nodes layer of load application in the z direction
        nzf = 1 + floor(zf/str.ew) + (mod(zf,str.ew) > str.ew/2); % Final nodes layer of load application in the z direction
        node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Initial node of load application
        
        if(nz0 == nzf)
            loads.nodes(p) = node; 
            loads.Fx(p) = (zf-z0)*Fx;
            loads.Fy(p) = (zf-z0)*Fy;
            loads.Fz(p) = (zf-z0)*Fz;
            p = p+1;
        else
            loads.nodes(p) = node; 
            loads.Fx(p) = ((nz0-1/2)*str.ew - z0)*Fx;
            loads.Fy(p) = ((nz0-1/2)*str.ew - z0)*Fy;
            loads.Fz(p) = ((nz0-1/2)*str.ew - z0)*Fz;
            p = p+1;
            node = node + str.nx*str.ny; 
            for i = (nz0+1):(nzf-1)
                loads.nodes(p) = node; 
                loads.Fx(p) = str.ew*Fx;
                loads.Fy(p) = str.ew*Fy;
                loads.Fz(p) = str.ew*Fz;
                p = p+1; 
                node = node + str.nx*str.ny;
            end
            loads.nodes(p) = node; 
            loads.Fx(p) = (zf - (nzf-3/2)*str.ew)*Fx;
            loads.Fy(p) = (zf - (nzf-3/2)*str.ew)*Fy;
            loads.Fz(p) = (zf - (nzf-3/2)*str.ew)*Fz;
            p = p+1;
        end
        
        str.loads = [str.loads; [x0,x0,y0,y0,z0,zf,Fx,Fy,Fz]]; % Save informations about the load
        continue;
    end 
    
    % Load distributed in an area parallel to the xy plane %
    if(opx == 1 && opy == 1 && opz == 0) 
        x0 = input('Initial x-coordinate of load application: ');
        xf = input('Final x-coordinate of load application: ');
        y0 = input('Initial y-coordinate of load application: ');
        yf = input('Final y-coordinate of load application: ');
        z0 = input('z-coordinate of load application: ');
        if(x0 < 0 || xf < x0 || y0 < 0 || yf < y0 || z0 < 0 || x0 > str.l || xf > str.l || y0 > str.h || yf > str.h || z0 > str.w) 
            fprintf('Error: the load is applied outside the domain.');
            continue;
        end
        
        nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Initial nodes layer of load application in the x direction 
        nxf = 1 + floor(xf/str.el) + (mod(xf,str.el) > str.el/2); % Final nodes layer of load application in the x direction
        ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Initial nodes layer of load application in the y direction
        nyf = 1 + floor(yf/str.eh) + (mod(yf,str.eh) > str.eh/2); % Final nodes layer of load application in the y direction
        nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Nodes layer of load application in the z direction
        node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Initial node of load application
        
        if(nx0 == nxf)
            d = xf-x0;
            if(ny0 == nyf)
                loads.nodes(p) = node;
                loads.Fx(p) = d*(yf-y0)*Fx; 
                loads.Fy(p) = d*(yf-y0)*Fy;
                loads.Fz(p) = d*(yf-y0)*Fz;
                p = p+1; 
            else
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((ny0-1/2)*str.eh - y0)*Fx;
                loads.Fy(p) = d*((ny0-1/2)*str.eh - y0)*Fy;
                loads.Fz(p) = d*((ny0-1/2)*str.eh - y0)*Fz;
                p = p+1;
                node = node + 1; 
                for i = (ny0+1):(nyf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.eh*Fx;
                    loads.Fy(p) = d*str.eh*Fy;
                    loads.Fz(p) = d*str.eh*Fz;
                    p = p+1; 
                    node = node + 1;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                loads.Fy(p) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                loads.Fz(p) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                p = p+1;
            end 
        else
            if(ny0 == nyf)
                d = yf-y0;
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((nx0-1/2)*str.el - x0)*Fx;
                loads.Fy(p) = d*((nx0-1/2)*str.el - x0)*Fy;
                loads.Fz(p) = d*((nx0-1/2)*str.el - x0)*Fz;
                p = p+1;
                node = node + str.ny; 
                for i = (nx0+1):(nxf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.el*Fx;
                    loads.Fy(p) = d*str.el*Fy;
                    loads.Fz(p) = d*str.el*Fz;
                    p = p+1; 
                    node = node + str.ny;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(xf - (nxf-3/2)*str.el)*Fx;
                loads.Fy(p) = d*(xf - (nxf-3/2)*str.el)*Fy;
                loads.Fz(p) = d*(xf - (nxf-3/2)*str.el)*Fz;
                p = p+1;
            else 
                d = (nx0-1/2)*str.el - x0;
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((ny0-1/2)*str.eh - y0)*Fx;
                loads.Fy(p) = d*((ny0-1/2)*str.eh - y0)*Fy;
                loads.Fz(p) = d*((ny0-1/2)*str.eh - y0)*Fz;
                p = p+1;
                node = node + 1; 
                for i = (ny0+1):(nyf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.eh*Fx;
                    loads.Fy(p) = d*str.eh*Fy;
                    loads.Fz(p) = d*str.eh*Fz;
                    p = p+1; 
                    node = node + 1;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                loads.Fy(p) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                loads.Fz(p) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                p = p+1;
                node = node + (str.ny - nyf) + ny0;
                
                d = str.el;
                for j = (nx0+1):(nxf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*((ny0-1/2)*str.eh - y0)*Fx;
                    loads.Fy(p) = d*((ny0-1/2)*str.eh - y0)*Fy;
                    loads.Fz(p) = d*((ny0-1/2)*str.eh - y0)*Fz;
                    p = p+1;
                    node = node + 1; 
                    for i = (ny0+1):(nyf-1)
                        loads.nodes(p) = node; 
                        loads.Fx(p) = d*str.eh*Fx;
                        loads.Fy(p) = d*str.eh*Fy;
                        loads.Fz(p) = d*str.eh*Fz;
                        p = p+1; 
                        node = node + 1;
                    end
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                    loads.Fy(p) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                    loads.Fz(p) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                    p = p+1;
                    node = node + (str.ny - nyf) + ny0;
                end
                
                d = xf - (nxf-3/2)*str.el;
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((ny0-1/2)*str.eh - y0)*Fx;
                loads.Fy(p) = d*((ny0-1/2)*str.eh - y0)*Fy;
                loads.Fz(p) = d*((ny0-1/2)*str.eh - y0)*Fz;
                p = p+1;
                node = node + 1; 
                for i = (ny0+1):(nyf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.eh*Fx;
                    loads.Fy(p) = d*str.eh*Fy;
                    loads.Fz(p) = d*str.eh*Fz;
                    p = p+1; 
                    node = node + 1;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                loads.Fy(p) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                loads.Fz(p) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                p = p+1;
            end 
        end
        
        str.loads = [str.loads; [x0,xf,y0,yf,z0,z0,Fx,Fy,Fz]]; % Save informations about the load
        continue;
    end 
    
    % Load distributed in an area parallel to the xz plane % 
    if(opx == 1 && opy == 0 && opz == 1) 
        x0 = input('Initial x-coordinate of load application: ');
        xf = input('Final x-coordinate of load application: ');
        y0 = input('y-coordinate of load application: ');
        z0 = input('Initial z-coordinate of load application: ');
        zf = input('Final z-coordinate of load application: ');

        if(x0 < 0 || xf < x0 || y0 < 0 || z0 < 0 || zf < z0 || x0 > str.l || xf > str.l || y0 > str.h || z0 > str.w || zf > str.w) 
            fprintf('Error: the load is applied outside the domain.');
            continue;
        end
        
        nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Initial nodes layer of load application in the x direction 
        nxf = 1 + floor(xf/str.el) + (mod(xf,str.el) > str.el/2); % Final nodes layer of load application in the x direction
        ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Nodes layer of load application in the y direction      
        nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Initial nodes layer of load application in the z direction
        nzf = 1 + floor(zf/str.ew) + (mod(zf,str.ew) > str.ew/2); % Final nodes layer of load application in the z direction
        node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Initial node of load application
        
        if(nx0 == nxf)
            d = xf-x0;
            if(nz0 == nzf)
                loads.nodes(p) = node;
                loads.Fx(p) = d*(zf-z0)*Fx; 
                loads.Fy(p) = d*(zf-z0)*Fy;
                loads.Fz(p) = d*(zf-z0)*Fz;
                p = p+1; 
            else
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((nz0-1/2)*str.ew - z0)*Fx;
                loads.Fy(p) = d*((nz0-1/2)*str.ew - z0)*Fy;
                loads.Fz(p) = d*((nz0-1/2)*str.ew - z0)*Fz;
                p = p+1;
                node = node + str.nx*str.ny; 
                for i = (nz0+1):(nzf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.ew*Fx;
                    loads.Fy(p) = d*str.ew*Fy;
                    loads.Fz(p) = d*str.ew*Fz;
                    p = p+1; 
                    node = node + str.nx*str.ny;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(zf - (nzf-3/2)*str.ew)*Fx;
                loads.Fy(p) = d*(zf - (nzf-3/2)*str.ew)*Fy;
                loads.Fz(p) = d*(zf - (nzf-3/2)*str.ew)*Fz;
                p = p+1;
            end 
        else
            if(nz0 == nzf)
                d = zf-z0;
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((nx0-1/2)*str.el - x0)*Fx;
                loads.Fy(p) = d*((nx0-1/2)*str.el - x0)*Fy;
                loads.Fz(p) = d*((nx0-1/2)*str.el - x0)*Fz;
                p = p+1;
                node = node + str.ny; 
                for i = (nx0+1):(nxf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.el*Fx;
                    loads.Fy(p) = d*str.el*Fy;
                    loads.Fz(p) = d*str.el*Fz;
                    p = p+1; 
                    node = node + str.ny;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(xf - (nxf-3/2)*str.el)*Fx;
                loads.Fy(p) = d*(xf - (nxf-3/2)*str.el)*Fy;
                loads.Fz(p) = d*(xf - (nxf-3/2)*str.el)*Fz;
                p = p+1;
            else 
                d = (nz0-1/2)*str.ew - z0;
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((nx0-1/2)*str.el - x0)*Fx;
                loads.Fy(p) = d*((nx0-1/2)*str.el - x0)*Fy;
                loads.Fz(p) = d*((nx0-1/2)*str.el - x0)*Fz;
                p = p+1;
                node = node + str.ny; 
                for i = (nx0+1):(nxf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.el*Fx;
                    loads.Fy(p) = d*str.el*Fy;
                    loads.Fz(p) = d*str.el*Fz;
                    p = p+1; 
                    node = node + str.ny;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(xf - (nxf-3/2)*str.el)*Fx;
                loads.Fy(p) = d*(xf - (nxf-3/2)*str.el)*Fy;
                loads.Fz(p) = d*(xf - (nxf-3/2)*str.el)*Fz;
                p = p+1;
                node = node + str.nx*str.ny - (nxf-nx0)*str.ny;
                
                d = str.ew;
                for j = (nz0+1):(nzf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*((nx0-1/2)*str.el - x0)*Fx;
                    loads.Fy(p) = d*((nx0-1/2)*str.el - x0)*Fy;
                    loads.Fz(p) = d*((nx0-1/2)*str.el - x0)*Fz;
                    p = p+1;
                    node = node + str.ny; 
                    for i = (nx0+1):(nxf-1)
                        loads.nodes(p) = node; 
                        loads.Fx(p) = d*str.el*Fx;
                        loads.Fy(p) = d*str.el*Fy;
                        loads.Fz(p) = d*str.el*Fz;
                        p = p+1; 
                        node = node + str.ny;
                    end
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*(xf - (nxf-3/2)*str.el)*Fx;
                    loads.Fy(p) = d*(xf - (nxf-3/2)*str.el)*Fy;
                    loads.Fz(p) = d*(xf - (nxf-3/2)*str.el)*Fz;
                    p = p+1;
                    node = node + str.nx*str.ny - (nxf-nx0)*str.ny;
                end
                
                d = zf - (nzf-3/2)*str.ew;
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((nx0-1/2)*str.el - x0)*Fx;
                loads.Fy(p) = d*((nx0-1/2)*str.el - x0)*Fy;
                loads.Fz(p) = d*((nx0-1/2)*str.el - x0)*Fz;
                p = p+1;
                node = node + str.ny; 
                for i = (nx0+1):(nxf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.el*Fx;
                    loads.Fy(p) = d*str.el*Fy;
                    loads.Fz(p) = d*str.el*Fz;
                    p = p+1; 
                    node = node + str.ny;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(xf - (nxf-3/2)*str.el)*Fx;
                loads.Fy(p) = d*(xf - (nxf-3/2)*str.el)*Fy;
                loads.Fz(p) = d*(xf - (nxf-3/2)*str.el)*Fz;
                p = p+1;
            end 
        end
        
        str.loads = [str.loads; [x0,xf,y0,y0,z0,zf,Fx,Fy,Fz]]; % Save informations about the load
        continue;
    end 
    
    % Load distributed in an area parallel to the yz plane %
    if(opx == 0 && opy == 1 && opz == 1) 
        x0 = input('x-coordinate of load application: ');
        y0 = input('Initial y-coordinate of load application: ');
        yf = input('Final y-coordinate of load application: ');
        z0 = input('Initial z-coordinate of load application: ');
        zf = input('Final z-coordinate of load application: ');
        if(x0 < 0 || y0 < 0 || yf < y0 || z0 < 0 || zf < z0 || x0 > str.l || y0 > str.h || yf > str.h || z0 > str.w || zf > str.w) 
            fprintf('Error: the load is applied outside the domain.');
            continue;
        end
        
        nx0 = 1 + floor(x0/str.el) + (mod(x0,str.el) > str.el/2); % Nodes layer of load application in the x direction 
        ny0 = 1 + floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2); % Initial nodes layer of load application in the y direction
        nyf = 1 + floor(yf/str.eh) + (mod(yf,str.eh) > str.eh/2); % Final nodes layer of load application in the y direction
        nz0 = 1 + floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2); % Initial nodes layer of load application in the z direction
        nzf = 1 + floor(zf/str.ew) + (mod(zf,str.ew) > str.ew/2); % Final nodes layer of load application in the z direction
        node = (nz0-1)*str.nx*str.ny + (nx0-1)*str.ny + ny0; % Initial node of load application
        
        if(nz0 == nzf)
            d = zf-z0;
            if(ny0 == nyf)
                loads.nodes(p) = node;
                loads.Fx(p) = d*(yf-y0)*Fx; 
                loads.Fy(p) = d*(yf-y0)*Fy;
                loads.Fz(p) = d*(yf-y0)*Fz;
                p = p+1; 
            else
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((ny0-1/2)*str.eh - y0)*Fx;
                loads.Fy(p) = d*((ny0-1/2)*str.eh - y0)*Fy;
                loads.Fz(p) = d*((ny0-1/2)*str.eh - y0)*Fz;
                p = p+1;
                node = node + 1; 
                for i = (ny0+1):(nyf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.eh*Fx;
                    loads.Fy(p) = d*str.eh*Fy;
                    loads.Fz(p) = d*str.eh*Fz;
                    p = p+1; 
                    node = node + 1;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                loads.Fy(p) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                loads.Fz(p) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                p = p+1;
            end 
        else
            if(ny0 == nyf)
                d = yf-y0;
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((nz0-1/2)*str.ew - z0)*Fx;
                loads.Fy(p) = d*((nz0-1/2)*str.ew - z0)*Fy;
                loads.Fz(p) = d*((nz0-1/2)*str.ew - z0)*Fz;
                p = p+1;
                node = node + str.nx*str.ny; 
                for i = (nz0+1):(nzf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.ew*Fx;
                    loads.Fy(p) = d*str.ew*Fy;
                    loads.Fz(p) = d*str.ew*Fz;
                    p = p+1; 
                    node = node + str.nx*str.ny;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(zf - (nzf-3/2)*str.ew)*Fx;
                loads.Fy(p) = d*(zf - (nzf-3/2)*str.ew)*Fy;
                loads.Fz(p) = d*(zf - (nzf-3/2)*str.ew)*Fz;
                p = p+1;
            else 
                d = (nz0-1/2)*str.ew - z0;
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((ny0-1/2)*str.eh - y0)*Fx;
                loads.Fy(p) = d*((ny0-1/2)*str.eh - y0)*Fy;
                loads.Fz(p) = d*((ny0-1/2)*str.eh - y0)*Fz;
                p = p+1;
                node = node + 1; 
                for i = (ny0+1):(nyf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.eh*Fx;
                    loads.Fy(p) = d*str.eh*Fy;
                    loads.Fz(p) = d*str.eh*Fz;
                    p = p+1; 
                    node = node + 1;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                loads.Fy(p) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                loads.Fz(p) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                p = p+1;
                node = node + str.nx*str.ny - (nyf-ny0);
                
                d = str.ew;
                for j = (nz0+1):(nzf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*((ny0-1/2)*str.eh - y0)*Fx;
                    loads.Fy(p) = d*((ny0-1/2)*str.eh - y0)*Fy;
                    loads.Fz(p) = d*((ny0-1/2)*str.eh - y0)*Fz;
                    p = p+1;
                    node = node + 1; 
                    for i = (ny0+1):(nyf-1)
                        loads.nodes(p) = node; 
                        loads.Fx(p) = d*str.eh*Fx;
                        loads.Fy(p) = d*str.eh*Fy;
                        loads.Fz(p) = d*str.eh*Fz;
                        p = p+1; 
                        node = node + 1;
                    end
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                    loads.Fy(p) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                    loads.Fz(p) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                    p = p+1;
                    node = node + str.nx*str.ny - (nyf-ny0);
                end
                
                d = zf - (nzf-3/2)*str.ew;
                loads.nodes(p) = node; 
                loads.Fx(p) = d*((ny0-1/2)*str.eh - y0)*Fx;
                loads.Fy(p) = d*((ny0-1/2)*str.eh - y0)*Fy;
                loads.Fz(p) = d*((ny0-1/2)*str.eh - y0)*Fz;
                p = p+1;
                node = node + 1; 
                for i = (ny0+1):(nyf-1)
                    loads.nodes(p) = node; 
                    loads.Fx(p) = d*str.eh*Fx;
                    loads.Fy(p) = d*str.eh*Fy;
                    loads.Fz(p) = d*str.eh*Fz;
                    p = p+1; 
                    node = node + 1;
                end
                loads.nodes(p) = node; 
                loads.Fx(p) = d*(yf - (nyf-3/2)*str.eh)*Fx;
                loads.Fy(p) = d*(yf - (nyf-3/2)*str.eh)*Fy;
                loads.Fz(p) = d*(yf - (nyf-3/2)*str.eh)*Fz;
                p = p+1;
            end 
        end
        
        str.loads = [str.loads; [x0,x0,y0,yf,z0,zf,Fx,Fy,Fz]]; % Save informations about the load
        continue;
    end
    
    % Volumetric load % 
    if(opx == 1 && opy == 1 && opz == 1)
        fprintf('Error: please, do not enter a volumetric load.'); 
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

% Void Regions (elements with density 0) % 
p = 1;
fixDens = zeros(str.nelem,1);
fixVal = zeros(str.nelem,1);
while(1)
    fprintf('\n');
    op = input('Continue adding void regions? (1 = Yes, 0 = No) ');
    if(op == 0)
        break;
    elseif(isempty(op) || op ~= 1)
        continue; 
    end
    
    x0 = input('Initial x-coordinate of the void region: '); 
    xf = input('Final x-coordinate of the void region: ');
    if((x0 < 0 && xf < 0) || (x0 > str.l && xf > str.l))
        fprintf('Error: void region is outside the domain.\n');
        continue;
    end
    y0 = input('Initial y-coordinate of the void region: '); 
    yf = input('Final y-coordinate of the void region: '); 
    if((y0 < 0 && yf < 0) || (y0 > str.h && yf > str.h))
        fprintf('Error: void region is outside the domain.\n');
        continue;
    end
    z0 = input('Initial z-coordinate of the void region: '); 
    zf = input('Final z-coordinate of the void region: '); 
    if((z0 < 0 && zf < 0) || (z0 > str.w && zf > str.w))
        fprintf('Error: void region is outside the domain.\n');
        continue;
    end
       
    ex0 = max(floor(x0/str.el) + (mod(x0,str.el) > str.el/2), 1); % Initial element layer of the fixed densities in the x direction 
    exf = min(floor(xf/str.el) + (mod(xf,str.el) > str.el/2), str.nelx); % Final element layer of the fixed densities in the x direction
    ey0 = max(floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2), 1); % Initial element layer of the fixed densities in the y direction
    eyf = min(floor(yf/str.eh) + (mod(yf,str.eh) > str.eh/2), str.nely); % Final element layer of the fixed densities in the y direction
    ez0 = max(floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2), 1); % Initial element layer of the fixed densities in the z direction
    ezf = min(floor(zf/str.ew) + (mod(zf,str.ew) > str.ew/2), str.nelz); % Final element layer of the fixed densities in the z direction    
    for elz = ez0:ezf
        for elx = ex0:exf
            for ely = ey0:eyf
                el = ely + (elx-1)*str.nely + (elz-1)*str.nelx*str.nely;
                fixDens(el) = 1;
            end
        end
    end 
    
    str.fix(p,:) = [ex0,exf,ey0,eyf,ez0,ezf,0]; 
    p = p + 1;
end

% Solid Regions (elements with density 1) % 
while(1)
    fprintf('\n');
    op = input('Continue adding solid regions? (1 = Yes, 0 = No) ');
    if(op == 0)
        break;
    elseif(isempty(op) || op ~= 1)
        continue; 
    end
    
    x0 = input('Initial x-coordinate of the solid region: '); 
    xf = input('Final x-coordinate of the solid region: ');
    if((x0 < 0 && xf < 0) || (x0 > str.l && xf > str.l))
        fprintf('Error: solid region is outside the domain.\n');
        continue;
    end
    y0 = input('Initial y-coordinate of the solid region: '); 
    yf = input('Final y-coordinate of the solid region: '); 
    if((y0 < 0 && yf < 0) || (y0 > str.h && yf > str.h))
        fprintf('Error: solid region is outside the domain.\n');
        continue;
    end
    z0 = input('Initial z-coordinate of the solid region: '); 
    zf = input('Final z-coordinate of the solid region: '); 
    if((z0 < 0 && zf < 0) || (z0 > str.w && zf > str.w))
        fprintf('Error: solid region is outside the domain.\n');
        continue;
    end
       
    ex0 = max(floor(x0/str.el) + (mod(x0,str.el) > str.el/2), 1); % Initial element layer of the fixed densities in the x direction 
    exf = min(floor(xf/str.el) + (mod(xf,str.el) > str.el/2), str.nelx); % Final element layer of the fixed densities in the x direction
    ey0 = max(floor(y0/str.eh) + (mod(y0,str.eh) > str.eh/2), 1); % Initial element layer of the fixed densities in the y direction
    eyf = min(floor(yf/str.eh) + (mod(yf,str.eh) > str.eh/2), str.nely); % Final element layer of the fixed densities in the y direction
    ez0 = max(floor(z0/str.ew) + (mod(z0,str.ew) > str.ew/2), 1); % Initial element layer of the fixed densities in the z direction
    ezf = min(floor(zf/str.ew) + (mod(zf,str.ew) > str.ew/2), str.nelz); % Final element layer of the fixed densities in the z direction    
    for elz = ez0:ezf
        for elx = ex0:exf
            for ely = ey0:eyf
                el = ely + (elx-1)*str.nely + (elz-1)*str.nelx*str.nely;
                fixDens(el) = 1;
                fixVal(el) = 1;
            end
        end
    end 
    
    str.fix(p,:) = [ex0,exf,ey0,eyf,ez0,ezf,1]; 
    p = p + 1;
end

dens = 1:str.nelem;
str.freeDens = dens(fixDens == 0); % elements with densities not fixed
str.fixedDens = dens(fixDens == 1); % elements with densities fixed 
str.fixedDensVal = fixVal(str.fixedDens); % density values of the elements with fixed densities
% ---------- %

end