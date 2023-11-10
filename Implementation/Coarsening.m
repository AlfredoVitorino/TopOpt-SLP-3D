function [C,F] = Coarsening(nnode,lambda,S,SC)
% Coarsening constructs the coarse grid for algebraic multigrid. %
% INPUT: nnode - number of algebraic nodes on the fine grid.
%        lambda - each component 'i' of this vector contains the number of nodes that are strongly influenced by node 'i'.
%        S - each row 'i' of this matrix contains the indexes of the nodes that have a strong connection with node 'i'.
%        SC - each row 'i' of this matrix contains the values of the strong connections of node 'i'.
% OUTPUT: C - indexes of the nodes on the coarse grid.
%         F - indexes of the nodes that are only on the fine grid.
% ---------- %

U = 1:nnode;
C = zeros(nnode,1); 
F = zeros(nnode,1); 
oldlambda = lambda; 
minS = min(nonzeros(S),[],2); 
maxS = max(S,[],2);
[maxlambda,i] = max(lambda); 
while(max(U) ~= 0 && maxlambda ~= 0) 
    C(i) = i; 
    U(i) = 0;  
    lambda(i) = 0;
    nS = S(i,1:oldlambda(i)); 
    v = nS(ismembc(nS,nonzeros(U(minS(i):maxS(i))))); 
    U(v) = 0; 
    lambda(v) = 0; 
    F(v) = v;
    W = ismembc(S(v,:),nonzeros(U(min(minS(v)):max(maxS(v))))); 
    for p = 1:length(v)
        j = v(p);
        nS = S(j,1:oldlambda(j));
        w = nS(W(p,:));
        lambda(w) = lambda(w) + 1; 
    end
    [maxlambda,i] = max(lambda);
end
F = nonzeros(F);
aux = zeros(length(F),1);
ind = 0; 
V = ismembc(S(F,:),nonzeros(C)); 
z = find(all(V==0,2));
if(any(z))
    for k = 1:length(z)
        i = F(z(k));    
        nS = S(i,1:oldlambda(i)); 
        v = nS(V(z(k),:)); 
        if(isempty(v))  
            [w,colS,colF] = intersect(nS,F);  
            if(~isempty(w))
                s = SC(i,colS); 
                [~,j] = max(s); 
                C(w(j)) = w(j); 
                ind = ind + 1; 
                aux(ind) = colF(j);  
                V = ismembc(S(F,:),nonzeros(C)); 
            end
        end
    end
end
C = nonzeros(C); 
F(aux(1:ind)) = [];

end