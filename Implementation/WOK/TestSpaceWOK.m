function [V] = TestSpaceWOK(B,n,nt,V0,tolr,tolq,maxiter,omega,nv,P,kv,dens,p,emin,str,elem,mr)
% TestSpaceWOK constructs the test space matrix solving the generalized eigenvalue problem, using the RQMCG method, when the global stiffness matrix is not explicitly generated. %
% INPUT: B - smoother iteration matrix. 
%        n - number of algebraic nodes. 
%        nt - number of test space vectors. 
%        V0 - initial guess for the eigenvectors. 
%        tolr - tolerance for the RQMCG method.
%        tolq - tolerance for the RQMCG method.
%        maxiter - maximum number of iterations for the RQMCG method.
%        omega - smoother relaxation parameter.
%        nv - grid level. 
%        P - cell containing the prolongation matrices.
%        kv - element stiffness matrix, converted into a column vector. 
%        dens - element densities vector (filtered).
%        p - penalty parameter of the SIMP model.
%        emin - Young's modulus of the void material
%        str - structure with the problem data.
%        elem - structure with element characteristics (deg - polinomial degree, type - element type (1 = lagrange, 2 = serendipity).
%        mr - structure that contains parameters used by the multiresolution method.
% OUTPUT: V - test space matrix. 
% ---------- %

V = zeros(n,nt);
Z = zeros(n,nt); 
for j = 1:nt
    U = V(:,1:(j-1));
    z = Z(:,1:(j-1))';
    y = V0(:,j); 
    y = y/norm(y);
    x = y - U*(U'*B*y); 
    %xA = A*x; 
    xA = AciProd(x,str,nv,P,kv,dens,p,emin,elem,mr);
    xB = B*x; 
    gamma = x'*xA; 
    eta = x'*xB; 
    q = gamma/eta;
    r = xA - q*xB;
    nr = norm(r)/norm(xA);
    nq = 1; 
    k = 1; 
    while(nr > tolr && nq > tolq && k <= maxiter) 
        g = -2*r/eta; 
        s = (1/omega)*(B\g);
        if(k == 1)
            db = s;
        else
            beta = (g'*s)/(oldg'*olds); 
            db = s + beta*d;
        end
        d = db - U*(z*db); 
        %dA = A*d;
        dA = AciProd(d,str,nv,P,kv,dens,p,emin,elem,mr);
        dB = B*d; 
        a1 = d'*xA; a2 = d'*dA; a3 = d'*xB; a4 = d'*dB;
        a = a2*a3 - a1*a4;
        b = a2*eta - a4*gamma; 
        c = a1*eta - a3*gamma; 
        if(b < 0)
            alpha = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        else
            alpha = -(2*c)/(b + sqrt(b^2 - 4*a*c));
        end
        x = x + alpha*d; 
        xA = xA + alpha*dA;
        xB = xB + alpha*dB; 
        gamma = gamma + 2*a1*alpha + a2*alpha^2; 
        eta = eta + 2*a3*alpha + a4*alpha^2; 
        oldq = q; 
        q = gamma/eta;
        oldg = g; 
        olds = s;       
        r = xA - q*xB;      
        nr = norm(r)/norm(xA);
        nq = abs(q-oldq)/q; 
        k = k + 1;           
    end
    V(:,j) = x/sqrt(eta); 
    Z(:,j) = B*V(:,j); 
end

end