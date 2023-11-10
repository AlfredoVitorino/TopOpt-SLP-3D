function [P,freedofsc] = Prolongation(nxf,nyf,nzf,nxc,nyc,nzc,nc,freedofsf)
% Prolongation obtains the prolongation matrix for geometric multigrid. %
% INPUT: nxf - number of fine nodes in the length(x) direction. 
%        nyf - number of fine nodes in the height(y) direction. 
%        nzf - number of fine nodes in the width(z) direction. 
%        nxc - number of coarse nodes in the length(x) direction. 
%        nyc - number of coarse nodes in the height(y) direction. 
%        nzc - number of coarse nodes in the width(z) direction. 
%        nc - total number of coarse nodes. 
%        freedofsf - indexes of the free (without support) degrees of freedom for the fine nodes. 
% OUTPUT: P - prolongation matrix. 
%         freedofsc - indexes of the free (without support) degrees of freedom for the coarse nodes. 
% ---------- %

w = [1,0.5,0.25,0.125]; 
maxind = 81*nc; 
iP = zeros(maxind,1); 
jP = zeros(maxind,1);
sP = zeros(maxind,1);
freedofsc = zeros(1,3*nc);
aux = zeros(1,3*nc);
ind = 0; 
% Loop in the coordinates (x,y,z) of the coarse nodes %
for zc = 1:nzc 
    for xc = 1:nxc 
        for yc = 1:nyc 
            col = (zc-1)*nxc*nyc + (xc-1)*nyc + yc; % index of the node on the coarse grid
            xf = 2*xc-1; yf = 2*yc-1; zf = 2*zc-1; % coordinates of the node on the fine grid
            % Loop in the neighbor nodes on the fine grid %
            for z = max(zf-1,1):min(zf+1,nzf)
                for x = max(xf-1,1):min(xf+1,nxf)
                    for y = max(yf-1,1):min(yf+1,nyf)
                        row = (z-1)*nxf*nyf + (x-1)*nyf + y; % index of the neighbor on the fine grid 
                        val = (xf-x)^2+(yf-y)^2+(zf-z)^2 + 1; 
                        ind = ind+1; iP(ind) = 3*row-2; jP(ind) = 3*col-2; sP(ind) = w(val); 
                        ind = ind+1; iP(ind) = 3*row-1; jP(ind) = 3*col-1; sP(ind) = w(val);
                        ind = ind+1; iP(ind) = 3*row; jP(ind) = 3*col; sP(ind) = w(val);                   
                        if(val == 1)
                            % Free (without support) degrees of freedom for the coarse nodes % 
                            freedofsc(3*col-2) = 3*col-2; 
                            freedofsc(3*col-1) = 3*col-1;
                            freedofsc(3*col) = 3*col;
                            aux(3*col-2) = 3*row-2; 
                            aux(3*col-1) = 3*row-1;
                            aux(3*col) = 3*row;
                        end
                    end
                end
            end
        end
    end
end
P = sparse(iP(1:ind),jP(1:ind),sP(1:ind)); 
freedofsc = freedofsc(ismember(aux,freedofsf)); 
P = P(freedofsf,freedofsc); 

end