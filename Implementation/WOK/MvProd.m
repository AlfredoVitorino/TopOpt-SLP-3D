function y = MvProd(x,B)
% MvProd applies the diagonal preconditioner, without constructing the global stiffness matrix. % 
% INPUT: x - function argument (vector that will be multiplied by the inverse of the diagonal).
%        B - inverse of the diagonal, as a vector.  
% OUTPUT: y - product of the inverse diagonal by the vector x. 
% ---------- %

y = B.*x;

end