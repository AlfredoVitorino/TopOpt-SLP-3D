function [xtotal] = TotalStr(str,xstar,sym,prjop)
% TotalStr generates the complete element densities vector of the optimal structure, based on the symmetries. % 
% INPUT: str - structure with the problem data. 
%        xstar - element densities vector (only the symmetric part). 
%        sym - structure with some symmetry options. 
%        prjop - indicates if the densities were rounded to 0 or 1 at the end of the SLP algorithm.
% OUTPUT: xtotal - complete element densities vector.
% ---------- %

x = xstar;
if(sym.xy)
    x2 = zeros(str.nelx*str.nely, str.nelz); 
    for i = 1:str.nelz
        x2(:,i) = x((str.nelz-i)*str.nelx*str.nely+1:(str.nelz-i+1)*str.nelx*str.nely);
    end
    x = [x2(:); x];
    str.nelz = str.nelz*2;
end
if(sym.yz)
    aux = reshape([x(:),zeros(rem(str.nelx*str.nely - rem(numel(x),str.nelx*str.nely),str.nelx*str.nely),1)],str.nelx*str.nely,[])'; 
    aux2 = zeros(size(aux,1),size(aux,2));
    i = 1;
    for k = str.nelx:-1:1
        for j = 1:str.nely
            aux2(:,i) = aux(:,(k-1)*str.nely+j);
            i = i+1;
        end
    end
    aux2 = [aux2,aux]';
    x2 = aux2(:);
    str.nelx = str.nelx*2;
    str.nx = str.nelx+1;
    x = x2;
end
xtotal = x;
   
str.nelem = str.nelem*(sym.xy+sym.yz+sym.xz)*2;

Display3D(str,xtotal,prjop);

end