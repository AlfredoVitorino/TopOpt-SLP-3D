function [xfil] = ApplyFilterWOK(x,str,rmin,w,opfilter,trans,mr,strDens,strDsgn)
% ApplyfilterWOK aplies the filter to the elements densities without using the matrix W. % 
% INPUT: x - vector with the original elements densities.
%        str - structure with the problem data.
%        rmin - filter radius.
%        w - vector with the sum of the weight factors for each finite element. 
%        opfilter - filter option (0 = no filter, 1 = mean density filter).
%        trans - boolean parameter that indicates whether w is to be multiplied to x or to xfil.
%        mr - structure that contains parameters used by the multiresolution method.
%        strDens - structure with density elements grid data for multiresolution. 
%        strDsgn - structure with design variable elements grid data for multiresolution.
% OUTPUT: xfil - vector with the filtered elements densities. 
% ---------- % 
 
if(opfilter == 0) % No filter
    xfil = x;
elseif(opfilter == 1 || opfilter == 2) % Mean density filter  
    if(mr.op)
        el = strDens.el;
        eh = strDens.eh;
        ew = strDens.ew;
        nelx = strDens.nelx;
        nely = strDens.nely;
        nelz = strDens.nelz; 
        if(mr.n ~= mr.d && trans)
            nelem = strDsgn.nelem; 
        else
            nelem = strDens.nelem;
        end
        if(mr.n == mr.d)
            fix = zeros(strDsgn.nelem,1);
            fix(strDsgn.fixedDens) = 1; % Indicates if the element has fixed density
        else
            fixDens = zeros(strDens.nelem,1);
            fixDens(strDens.fixedDens) = 1; 
            fixDsgn = zeros(strDsgn.nelem,1);
            fixDsgn(strDsgn.fixedDens) = 1;
        end
    else
        el = str.el;
        eh = str.eh;
        ew = str.ew;
        nelx = str.nelx;
        nely = str.nely;
        nelz = str.nelz; 
        nelem = str.nelem; 
        fix = zeros(str.nelem,1);
        fix(str.fixedDens) = 1; % Indicates if the element has fixed density
    end
  
    xfil = zeros(nelem,1);
    
    if(trans)
        x = (x./w);
    end
    
    if(~mr.op || mr.n == mr.d)
        frx = floor(rmin/el);
        fry = floor(rmin/eh);
        frz = floor(rmin/ew);
        rmin2 = rmin^2;  
        relx = zeros(frz+1);
        rely = zeros(frz+1,frx+1);
        d2 = zeros(frx+1,fry+1,frz+1);
        if(opfilter == 1)
            alpha = 1/(2*((rmin/3)^2));
            beta = 1/(2*pi*(rmin/3));
        end

        for k1 = 1:frz+1
            dz2 = ((k1-1)*ew)^2;
            dx = floor(sqrt(rmin2-dz2)/el);
            relx(k1) = dx;
            for i1 = 1:dx+1
                dx2 = ((i1-1)*el)^2;
                dy = floor(sqrt(rmin2-dz2-dx2)/eh);
                rely(k1,i1)=dy;
                if(opfilter == 1)
                    d2(i1,1:dy+1,k1) = exp(-(dz2+dx2+((0:dy)*eh).^2)*alpha)*beta;
                elseif(opfilter == 2)
                    d2(i1,1:dy+1,k1) = (rmin-sqrt(dz2+dx2+((0:dy)*str.eh).^2))/rmin; 
                end
            end
        end

        % Multiplying the weight factors by x %
        e1 = 0;
        for k1 = 1:nelz
            for i1 = 1:nelx
                for j1 = 1:nely
                    e1 = e1 + 1;
                    if(fix(e1) == 0)
                        for k2 = k1:min(k1+frz, nelz)
                            pz = k2-k1+1;
                            e21 = (k2-1)*nelx*nely;
                            for i2 = max(i1-relx(pz), 1):min(i1+relx(pz), nelx)
                                px = abs(i2-i1)+1;
                                e22 = e21 + (i2-1)*nely;
                                minj2 = max(j1-rely(pz,px), 1);
                                maxj2 = min(j1+rely(pz,px), nely);
                                if (minj2 <= (e1-e22))
                                    minj2 = e1-e22+1;
                                    if (maxj2 >= (e1-e22))
                                        if(opfilter == 1)
                                            xfil(e1) = xfil(e1)+beta*x(e1);
                                        elseif(opfilter == 2)
                                            xfil(e1) = xfil(e1)+x(e1);
                                        end
                                    end
                                end
                                for j2 = minj2:maxj2
                                    py = abs(j2-j1)+1;
                                    e2 = e22 + j2;
                                    if(fix(e2) == 0)
                                        xfil(e1) = xfil(e1)+d2(px,py,pz)*x(e2);
                                        xfil(e2) = xfil(e2)+d2(px,py,pz)*x(e1);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    else % Multiresolution with n~=d
        frx = floor(rmin/strDens.el); 
        fry = floor(rmin/strDens.eh); 
        frz = floor(rmin/strDens.ew);
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
                    nbd = nb(dind);
                    for j3 = 1:length(d)
                        if(fixDens(e1) == 0 && fixDsgn(nbd(j3)) == 0 || fixDens(e1) == 1 && fixDsgn(nbd(j3)) == 1)
                            if(opfilter == 1)
                                if(trans)
                                    xfil(nbd(j3)) = xfil(nbd(j3)) + (exp(-(d(j3)^2)*alpha)*beta)*x(e1);
                                else
                                    xfil(e1) = xfil(e1) + (exp(-(d(j3)^2)*alpha)*beta)*x(nbd(j3));
                                end
                            elseif(opfilter == 2)
                                if(trans)
                                    xfil(nbd(j3)) = xfil(nbd(j3)) + ((rmin-d)/rmin)*x(e1);
                                else
                                    xfil(e1) = xfil(e1) + ((rmin-d)/rmin)*x(nbd(j3));
                                end
                            end
                        end
                    end
                    ind = ld; 
                end
            end
        end
    end
    
    if (~trans)
        xfil = xfil./w;
    end      
end

end