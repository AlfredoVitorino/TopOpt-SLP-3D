function P = ProlongationDPLS(nnode,S,V,C,F,itp,kappa)
% ProlongationDPLS obtains the prolongation matrix for algebraic multigrid, using the DPLS method. %
% INPUT: nnode - number of algebraic nodes on the fine grid.
%        S - each row 'i' of this matrix contains the indexes of the nodes that have a strong connection with node 'i'.
%        V - test space matrix. 
%        C - indexes of the nodes on the coarse grid.
%        F - indexes of the nodes that are only on the fine grid.
%        itp - maximum number of iterations on DPLS method. 
%        kappa - tolerance for the condition number on DPLS method.
% OUTPUT: P - prolongation matrix. 
% ---------- %

iP = zeros(81*nnode,1); 
jP = zeros(81*nnode,1); 
sP = zeros(81*nnode,1); 
ind = 0; 
for j = 1:length(C)
    ind = ind + 1; 
    iP(ind) = C(j); 
    jP(ind) = j; 
    sP(ind) = 1; 
end
[Lj,Lc] = ismember(S,C);
for j = 1:length(F)
    i = F(j); 
    J = S(i,Lj(i,:));
    if(~isempty(J))
        col = nonzeros(Lc(i,:));
        Vb = V(J,:)'; 
        r = V(i,:)'; 
        lj = length(J);
        R = []; 
        oldind = ind + 1;
        k = 0;
        while(k < lj && k < itp)
            k = k+1; 
            s = ((Vb'*r).^2)./((r'*r)*sum(Vb.*Vb)');
            [~,maxind] = max(s);
            oldR = R; 
            [Q,R] = qr([R, Vb(:,maxind)]);
            cd = 1/rcond(R(1:k,1:k));
            if(cd <= kappa)
                J(maxind) = [];
                Vb(:,maxind) = [];
                r = Q*r;       
                Vb = Q*Vb; 
                ind = ind + 1; 
                iP(ind) = i;  
                jP(ind) = col(maxind); 
                col(maxind) = []; 
            else
                R = oldR; 
                break;
            end            
        end
        w = R\r; 
        sP(oldind:ind) = w';
    end
end
P = sparse(iP(1:ind),jP(1:ind),sP(1:ind),nnode,length(C));
end