function [V] = TestSpace(A,B,n,nt,V0,tolr,tolq,maxiter,opcond)
% TestSpace constructs the test space matrix solving the generalized eigenvalue problem, using the RQMCG method. %
% INPUT: A - system matrix. 
%        B - smoother iteration matrix. 
%        n - number of algebraic nodes. 
%        nt - number of test space vectors. 
%        V0 - initial guess for the eigenvectors. 
%        tolr - tolerance for the RQMCG method.
%        tolq - tolerance for the RQMCG method.
%        maxiter - maximum number of iterations for the RQMCG method.
%        opcond - preconditioner option (0 - diag, 1 - ichol).
% OUTPUT: V - test space matrix. 
% ---------- %

if(opcond == 0) % diag
    M = diag(diag(A));   
elseif(opcond == 1) % ichol
    shift = 0.1;
    M = ichol(A, struct('diagcomp', shift)); 
end

V = zeros(n,nt);
Z = zeros(n,nt); 
for j = 1:nt
    U = V(:,1:(j-1));
    z = Z(:,1:(j-1))';
    y = V0(:,j); 
    y = y/norm(y);
    x = y - U*(U'*B*y); 
    xA = A*x; 
    xB = B*x; 
    gamma = x'*xA; 
    eta = x'*xB; 
    q = gamma/eta;
    r = xA - q*xB;
    nr = norm(r)/norm(xA);
    nq = 1; 
    k = 1; 
    %fprintf('\nAutovetor %d\n', j); 
    %fprintf('\nit: %4d | nr: %d | nq: %d \n', k, nr, nq);
    while(nr > tolr && nq > tolq && k <= maxiter) 
        g = -2*r/eta; 
        if(opcond == 0) %diag
            s = M\g;
        elseif(opcond == 1) %ichol
            s = M'\(M\g); 
        end
        if(k == 1)
            db = s;
        else
            beta = (g'*s)/(oldg'*olds); 
            db = s + beta*d;
        end
        d = db - U*(z*db); 
        dA = A*d; 
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
        %fprintf('it: %4d | nr: %d | nq: %d \n', k, nr, nq); 
    end
    V(:,j) = x/sqrt(eta); 
    Z(:,j) = B*V(:,j); 
end

end