function [iK,jK,perm] = SetupStiffMatrix(str,x,p,elem,Dofs,mr,strDens)
% SetupStiffMatrix obtains the row and column indexes with nonzero elements of the global stiffness matrix and its permutation vector. % 
% INPUT: str - structure with the problem data.  
%        x - initial element densities vector (filtered).
%        p - penalty parameter of the SIMP model.
%        elem - structure with element characteristics (deg - polynomial degree, type - element type (1 = Lagrange, 2 = serendipity).
%        Dofs - matrix with the indexes of the degrees of freedom for each element.
%        mr - structure that contains parameters used by the multiresolution method.
%        strDens - structure with density mesh data for multiresolution. 
% OUTPUT: iK - row indexes with nonzero elements of the matrix K.
%         jK - column indexes with nonzero elements of the matrix K.
%         perm - approximate minimum degree permutation vector of K.
% ---------- % 

% m = total number of nodes in each element
if(elem.type == 1) % Lagrange element
    m = (elem.deg+1)^3;
elseif(elem.type == 2) % Serendipity element
    m = 8 + 12*(elem.deg-1);
end

iK = kron(Dofs, ones(3*m,1))';
iK = iK(:);
jK = kron(Dofs, ones(1,3*m))';
jK = jK(:);

if(~mr.op)
    sK = (100*rand((3*m)^2,1))*(x'.^p); 
    sK = sK(:);
else % Multiresolution
    sK = zeros(size(iK,1),1); 
    indS = 1; 
    j = 1; 
    nely = strDens.nely; 
    nelxy = strDens.nely*strDens.nelx;
    indstep = (3*m)^2;
    for i = 1:str.nelem   
        vb = zeros(mr.n^2, 1); 
        ind = 1; 
        for n1 = 0:(mr.n-1) 
            for n2 = 0:(mr.n-1)
                vb(ind) = j + n1*nelxy + n2*nely; 
                ind = ind + 1; 
            end
        end     
        v = zeros(mr.n^2, mr.n);
        for n3 = 0:(mr.n-1)
            v(:,n3+1) = vb+n3; 
        end
        v = v(:);         
        v = sort(v);

        xel = x(v).^p;
        sK(indS:(indS+indstep-1)) = sum(xel)*(100*rand(indstep,1)); 
        indS = indS+indstep; 

        j = j+mr.n;
        if (mod(j-1,nely) == 0)
            j = j+(mr.n-1)*nely;
            if (mod(j,nelxy) == 1)
                j = j+(mr.n-1)*nelxy;
            end
        end
    end
end

K = sparse(iK,jK,sK);
K = K(str.freedofs,str.freedofs);
perm = amd(K); 

end