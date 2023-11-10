function [Fgrad] = Grad(str,u,p,emin,k,x,Dofs,W,w,opfilter,mr,strDens)
% Grad calculates the gradient vector of the objective function of the topology optimization problem. %
% INPUT: str - structure with the problem data. 
%        u - nodal displacements vector. 
%        p - penalty parameter for the SIMP model. 
%        emin - Young's modulus of the void material.
%        k - finite element stiffness matrix. 
%        x - element densities vector (filtered). 
%        Dofs - matrix with the indexes of the degrees of freedom for each element. 
%        W - each row of this matrix contains the weight factors for each finite element, associated to the average density filter. 
%        w - vector with the sum of the weight factors for each finite element. 
%        opfilter - filter option (0 = no filter, 1 = weighted average density filter, 2 = average density filter).
%        mr - structure that contains parameters used by the multiresolution method.
%        strDens - structure with density mesh data for multiresolution.
% OUTPUT: Fgrad - gradient vector of the objective function. 
% ---------- % 

if(~mr.op)
    Fgrad0 = zeros(str.nelem,1); 
    for i = str.freeDensG
        uel = u(Dofs(i,:));
        Fgrad0(i) = -p*(x(i)^(p-1))*(str.E-emin)*(uel'*k*uel);        
    end
else % Multiresolution
    Fgrad0 = zeros(strDens.nelem,1);
    nely = strDens.nely; 
    nelxy = strDens.nelx*strDens.nely;
    for i = str.freeDensG
        vb = zeros(mr.n^2, 1); 
        ind = 1; 
        for n1 = 0:(mr.n-1) 
            for n2 = 0:(mr.n-1)
                vb(ind) = strDens.fde(i) + n1*nelxy + n2*nely; 
                ind = ind + 1; 
            end
        end     
        v = zeros(mr.n^2, mr.n);
        for n3 = 0:(mr.n-1)
            v(:,n3+1) = vb+n3; 
        end
        v = v(:);         
        v = sort(v);
        
        uel = u(Dofs(i,:)); % Nodal displacements of the finite element i 
        if (~mr.interp) % Using the original displacements to compute the gradient. 
            for i2 = 1:length(v)  
                kel = k{i2}; % Stiffness integrand for the density element
                Fgrad0(v(i2)) = -p*(x(v(i2))^(p-1))*(str.E-emin)*(uel'*kel*uel);             
            end
        else % Using interpolated displacements to compute the gradient.
            ud = mr.Idisp*uel; % ud contains the interpolated displacements of all density elements inside the finite element i
            for i2 = 1:length(v) 
                ueli = ud(mr.Ldofs(i2,:)); % Interpolated displacements of a single density element         
                %kel = k{i2}; % Stiffness integrand for the density element
                kel = mr.kd;
                Fgrad0(v(i2)) = -p*(x(v(i2))^(p-1))*(str.E-emin)*(ueli'*kel*ueli);             
            end
        end
    end
end

Fgrad = ApplyFilter(Fgrad0,W,w,opfilter,true);

end