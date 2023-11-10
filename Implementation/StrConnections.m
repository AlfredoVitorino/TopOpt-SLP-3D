function [S,SC,lambda] = StrConnections(A,V,nnode,theta,optheta)
% StrConnections obtains the strength of connections for each algebraic node. %
% INPUT: A - system matrix. 
%        V - test space matrix. 
%        nnode - number of algebraic nodes. 
%        theta - strong connections parameter.
%        optheta - option if 'theta' will be used as (0 - a tolerance, 1 - the maximum number of strong connections per node).
% OUTPUT: S - each row 'i' of this matrix contains the indexes of the nodes that have a strong connection with node 'i'.
%         SC - each row 'i' of this matrix contains the values of the strong connections of node 'i'.
%         lambda - each component 'i' of this vector contains the number of nodes that are strongly influenced by node 'i'.
% ---------- %

S = zeros(nnode,80); 
SC = zeros(nnode,80); 
lambda = zeros(nnode,1); 
[iB,jB,~] = find(A - diag(diag(A)));
[~,uB] = unique(jB); 
indB = [uB; length(jB)+1];  
alpha = sum(V.*V,2);
for i = 1:nnode
    c = iB(indB(i):(indB(i+1)-1)); 
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