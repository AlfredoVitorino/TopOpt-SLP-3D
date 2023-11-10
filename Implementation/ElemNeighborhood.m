function [W,w] = ElemNeighborhood(str,rmin,opfilter)
% ElemNeighborhood finds the neighborhood and the weight factors of each finite element, for the filter application. % 
% INPUT: str - structure with the problem data. 
%        rmin - filter radius. 
%        opfilter - filter option (0 = no filter, 1 = weighted average density filter, 2 = average density filter).
% OUTPUT: W - each row of this matrix contains the weight factors for each finite element, associated to the average density filter. 
%         w - vector with the sum of the weight factors for each finite element. 
% ---------- %

frx = floor(rmin/str.el);
fry = floor(rmin/str.eh);
frz = floor(rmin/str.ew);
rmin2 = rmin^2;
if(opfilter == 1)
    alpha = 1/(2*((rmin/3)^2));
    beta = 1/(2*pi*(rmin/3));
end
iW = zeros(str.nelem,1); 
jW = zeros(str.nelem,1); 
sW = zeros(str.nelem,1); 
relx = zeros(frz+1);
rely = zeros(frz+1,frx+1);
d2 = zeros(frx+1,fry+1,frz+1);

for k1 = 1:frz+1
    dz2 = ((k1-1)*str.ew)^2;
    dx = floor(sqrt(rmin2-dz2)/str.el);
    relx(k1) = dx;
    for i1 = 1:dx+1
        dx2 = ((i1-1)*str.el)^2;
        dy = floor(sqrt(rmin2-dz2-dx2)/str.eh);
        rely(k1,i1) = dy;
        if(opfilter == 1)
            d2(i1,1:dy+1,k1) = exp(-(dz2+dx2+((0:dy)*str.eh).^2)*alpha)*beta;
        elseif(opfilter == 2)
            d2(i1,1:dy+1,k1) = (rmin-sqrt(dz2+dx2+((0:dy)*str.eh).^2))/rmin; 
        end
    end
end

e1 = 0;
ind = 1; 
for k1 = 1:str.nelz
    for i1 = 1:str.nelx
        for j1 = 1:str.nely
            e1 = e1 + 1;
            for k2 = k1:min(k1+frz, str.nelz)
                pz = k2-k1+1;
                e21 = (k2-1)*str.nelx*str.nely;
                for i2 = max(i1-relx(pz), 1):min(i1+relx(pz), str.nelx)
                    px = abs(i2-i1)+1;
                    e22 = e21 + (i2-1)*str.nely;
                    minj2 = max(j1-rely(pz,px), 1);
                    maxj2 = min(j1+rely(pz,px), str.nely);
                    if (minj2<=(e1-e22))
                        minj2 = e1-e22+1;
                        if (maxj2>=(e1-e22))
                            iW(ind) = e1; 
                            jW(ind) = e1; 
                            if(opfilter == 1)
                                sW(ind) = beta;
                            elseif(opfilter == 2)
                                sW(ind) = 1;
                            end
                            ind = ind + 1;
                        end
                    end
                    j2 = minj2:maxj2;
                    py = abs(j2-j1)+1;
                    e2 = e22 + j2;
                    wt = d2(px,py,pz);
                    indf = ind+length(wt)-1;
                    iW(ind:indf) = e1;
                    jW(ind:indf) = e2;
                    sW(ind:indf) = wt;
                    ind = indf + 1;
                    indf = ind+length(wt)-1;
                    iW(ind:indf) = e2;
                    jW(ind:indf) = e1;
                    sW(ind:indf) = wt;
                    ind = indf + 1;
                end
            end
        end
    end
end

W = sparse(iW,jW,sW);

if(~isempty(str.fixedDens))
    W(str.fixedDens,:) = 0.0; 
    W(:,str.fixedDens) = 0.0;
    newdiag = diag(W); 
    newdiag(str.fixedDens) = 1.0;
    W = W - diag(diag(W)) + diag(newdiag);
end

w = full(sum(W,2)); 

end