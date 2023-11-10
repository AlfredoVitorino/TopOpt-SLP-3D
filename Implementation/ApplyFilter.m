function [xfil] = ApplyFilter(x,W,w,opfilter,trans)
% Applyfilter applies the filter to the density elements. % 
% INPUT: x - vector with the original density elements.
%        W - each row of this matrix contains the weight factors for each element, associated to the average density filter. 
%        w - vector with the sum of the weight factors for each element. 
%        opfilter - filter option (0 = no filter, 1 = weighted average density filter, 2 = average density filter).
%        trans - boolean parameter that indicates whether w is to be multiplied to x or to xfil.
% OUTPUT: xfil - vector with the filtered density elements. 
% ---------- % 
 
if(opfilter == 0) % No filter
    xfil = x;
elseif(opfilter == 1 || opfilter == 2) % Average density filter
    if (trans)
        xfil = W'*(x./w);
    else
        xfil = (W*x)./w;
    end
end

end