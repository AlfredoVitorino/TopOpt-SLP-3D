function f = HeavRoot(eta,beta,x,velem,s)
% HeavRoot is a function used to find the parameter eta of the Heaviside projector. %
% INPUT: eta - current value for eta.
%        beta - second parameter of the Heaviside function.
%        x - element densities.
%        velem - element volumes.
%        s - current value of velem'*x.
% OUTPUT: f - value of the Heaviside function. 
% ---------- % 

c = tanh(beta*eta);
d = 1.0/(c+tanh(beta*(1.0-eta)));
f = velem'*((c+tanh(beta*(x-eta)))*d)-s;

end
