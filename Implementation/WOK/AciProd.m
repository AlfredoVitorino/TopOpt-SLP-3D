function y = AciProd(x,str,ind,P,kv,dens,p,emin,elem,mr)
% AciProd calculates the product Ac{ind}*x without constructing the matrix K. %
% INPUT: x - function argument (vector that will be multiplied by K).
%        str - structure with the problem data. 
%        ind - index of the Ac matrix.
%        P - prolongation matrix.
%        kv - element stiffness matrix, converted into a column vector.
%        dens - element densities vector (filtered). 
%        p - penalty parameter for the SIMP model. 
%        emin - Young's modulus of the void material.
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: y - product K*x.
% ---------- %

for i = ind-1:-1:1
    x = P{i}*x;  
end

y = KelemProd(str,x,kv,dens,p,emin,elem,mr);

for i = 1:ind-1
    y = P{i}'*y;  
end

end