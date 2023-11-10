function [W,w] = ElemNeighborhoodMR(strDens,strDsgn,rmin,mr,opfilter)
% ElemNeighborhoodMR finds the neighborhood between the density elements and design variable elements, for multiresolution. % 
% INPUT: strDens - structure with density mesh data for multiresolution. 
%        strDsgn - structure with design variable mesh data for multiresolution.
%        rmin - filter radius. 
%        mr - structure that contains parameters used by the multiresolution method.
%        opfilter - filter option (0 = no filter, 1 = weighted average density filter, 2 = average density filter).
% OUTPUT: W - each row of this matrix contains the weight factors for each density element.
%         w - vector with the sum of the weight factors for each density element. 
% ---------- %

frx = floor(rmin/strDens.el); 
fry = floor(rmin/strDens.eh); 
frz = floor(rmin/strDens.ew);
iW = zeros(strDens.nelem,1); 
jW = zeros(strDens.nelem,1); 
sW = zeros(strDens.nelem,1);
if(opfilter == 1)
    alpha = 1/(2*((rmin/3)^2)); 
    beta = 1/(2*pi*(rmin/3));
end
rmin2 = rmin^2;
elhDens = strDens.el/2;
ehhDens = strDens.eh/2;
ewhDens = strDens.ew/2;
elhDsgn = strDsgn.el/2;
ehhDsgn = strDsgn.eh/2;
ewhDsgn = strDsgn.ew/2;
nelxyDsgn = strDsgn.nelx*strDsgn.nely;
gridratio = mr.n/mr.d;
ind = 1;
e1 = 0;
for k1 = 1:strDens.nelz
    z1 = ewhDens + (k1-1)*strDens.ew;
    mink2 = max(round(max(k1-frz, 1)/gridratio), 1);
    maxk2 = min(round(min(k1+frz, strDens.nelz)/gridratio), strDsgn.nelz);
    k2 = mink2:maxk2;
    z2 = ewhDsgn + (k2-1)*strDsgn.ew;
    dz2 = (z2-z1).^2;
    for i1 = 1:strDens.nelx
        x1 = elhDens + (i1-1)*strDens.el;
        mini2 = max(round(max(i1-frx, 1)/gridratio), 1);
        maxi2 = min(round(min(i1+frx, strDens.nelx)/gridratio), strDsgn.nelx);
        i2 = mini2:maxi2;
        x2 = elhDsgn + (i2-1)*strDsgn.el;
        dx2 = (x2-x1).^2;
        for j1 = 1:strDens.nely
            e1 = e1 + 1;
            y1 = ehhDens + (j1-1)*strDens.eh;
            minj2 = max(round(max(j1-fry, 1)/gridratio), 1);
            maxj2 = min(round(min(j1+fry, strDens.nely)/gridratio), strDsgn.nely);   
            j2 = minj2:maxj2;
            y2 = ehhDsgn + (j2-1)*strDsgn.eh; 
            dy2 = (y2-y1).^2;
            d2 = repelem(dx2,length(dy2)) + repmat(dy2,1,length(dx2));
            d2 = repelem(dz2,length(d2)) + repmat(d2,1,length(dz2));            
            e21 = (k2-1)*nelxyDsgn;
            e22 = (i2-1)*strDsgn.nely;
            nb = repelem(e21,length(e22)) + repmat(e22,1,length(e21));
            nb = repelem(nb,length(j2)) + repmat(j2,1,length(nb)); 
            [~,dind] = find(d2 <= rmin2);
            d = sqrt(d2(dind));
            ld = ind + length(d);
            iW(ind:ld-1) = e1;
            jW(ind:ld-1) = nb(dind);
            if(opfilter == 1)
                sW(ind:ld-1) = exp(-(d.^2)*alpha)*beta;
            elseif(opfilter == 2)
                sW(ind:ld-1) = (rmin-d)/rmin;
            end
            ind = ld;          
        end
    end
end

W = sparse(iW,jW,sW); 

if(~isempty(strDsgn.fixedDens))
    W(strDens.fixedDens,strDsgn.freeDens) = 0.0; 
    W(strDens.freeDens,strDsgn.fixedDens) = 0.0;
end

w = full(sum(W,2)); 

end