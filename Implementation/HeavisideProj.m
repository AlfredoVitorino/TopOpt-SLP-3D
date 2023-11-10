function [xstar,prj] = HeavisideProj(x,prj,velem,scur,sdes,betaupd)
% HeavisideProj updates the beta and eta parameters used to compute the Heaviside projection of the density vector. %
% INPUT: x - vector with the original element densities.
%        prj - structure that contains the parameters of the Heaviside projector.
%        velem - element volumes.
%        scur - current structure volume.
%        sdes - desired structure volume.
%        s - desired structure volume.
%        betaupd - parameter that indicates if beta is to be updated or not.
% OUTPUT: xstar - vector with the projected densities.
%         prj - structure with updated eta and beta.
% ---------- % 

% Updating beta.
if (betaupd)
  prj.beta = min(prj.beta*prj.dbeta,prj.betamax);
end

% Ensuring that the values of f at eta = 0 and eta = 1 have opposite signs.
s = sdes;
f0 = HeavRoot(0.0,prj.beta,x,velem,s);
f1 = HeavRoot(1.0,prj.beta,x,velem,s);
i = 0;
while (((f0*f1)>0)&&(i<20))
   i = i+1;
   s = scur*(i/20.0)+sdes*(1-i/20.0);
  f0 = HeavRoot(0.0,prj.beta,x,velem,s);
  f1 = HeavRoot(1.0,prj.beta,x,velem,s);
end

% Computing eta.
eta = prj.eta;
fun = @(eta) HeavRoot(eta,prj.beta,x,velem,s);
prj.eta = fzero(fun,[0.0 1.0]);

% Updating x.
c = tanh(prj.beta*prj.eta);
cd = 1.0/(c+tanh(prj.beta*(1.0-prj.eta)));
xstar = (c+tanh(prj.beta*(x-prj.eta)))*cd;

end