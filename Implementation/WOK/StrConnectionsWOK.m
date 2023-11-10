function [S,SC,lambda] = StrConnectionsWOK(V,nnode,theta,optheta,nv,P,kv,dens,p,emin,str,elem,mr)
% StrConnectionsWOK obtains the strong connections for each algebraic node when the global stiffness matrix is not explicitly generated. %
% INPUT: V - test space matrix. 
%        nnode - number of algebraic nodes. 
%        theta - strong connections parameter.
%        optheta - option if 'theta' will be used as (0 - a tolerance, 1 - the maximum number of strong connections per node).
%        nv - grid level. 
%        P -  cell containing the prolongation matrices. 
%        kv - element stiffness matrix, converted into a column vector. 
%        dens - element densities vector (filtered).
%        p - penalty parameter of the SIMP model.
%        emin - Young's modulus of the void material
%        str - structure with the problem data. 
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: S - each row 'i' of this matrix contains the indexes of the nodes that have a strong connection with node 'i'.
%         SC - each row 'i' of this matrix contains the values of the strong connections of node 'i'.
%         lambda - each component 'i' of this vector contains the number of nodes that are strongly influenced by node 'i'.
% ---------- %

S = zeros(nnode,80); 
SC = zeros(nnode,80); 
lambda = zeros(nnode,1); 
alpha = sum(V.*V,2);
for i = 1:nnode
    e = zeros(nnode,1); 
    e(i) = 1; 
    a = AciProd(e,str,nv,P,kv,dens,p,emin,elem,mr);
    a(i) = 0;
    [c,~,~] = find(a);
    v = V(i,:)';
    W = V(c,:);
    z = W*v; 
    beta = z.*z;  
    gamma = alpha(i)*alpha(c);
    s = beta./gamma;
    if(optheta == 0) % Using theta as a tolerance 
        st = find(s > theta);
    else % Using theta as the maximum number of strong connections per node
        [~,spos] = sort(s,'descend'); 
        st = spos(1:min(theta,length(spos)));
    end
    ind = lambda(i)+1;
    lambda(i) = lambda(i) + length(st);
    S(i,ind:lambda(i)) = c(st);
    SC(i,ind:lambda(i)) = s(st);
end 

end