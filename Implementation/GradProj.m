function xnew = GradProj(x,g,prj,Vmax,velem)
% GradProj rounds variables to 0 or 1 according to the gradient of the Lagrangian. %
% INPUT: x - current approximate solution.
%        g - gradient of the Lagrangian.
%        prj - structure that contains the parameters used for rounding the optimal densities to 0 or 1
%        Vmax - upper limit for the volume of the structure.
%        velem - volume of each element of the domain.
% OUTPUT: xnew - new approximate solution. 
% ---------- %

g = -g; % Now, g stores -grad(Lag)
gnorm = norm(g,2);
coslim = cosd(prj.maxang); % Maximum value for the cosine of the angle between d and g
gt = prj.gthresh*norm(g,inf); % Gradient threshold. g(i) is considered to be near 0 if abs(g(i)) <= g
volfrac = (velem'*x)/sum(velem); % Fraction of the domain that is filled by the current solution x

xtop = find(x>=prj.ut); % indexes of the elements with density >= ut
ntop = length(xtop); 
done = round(Vmax/velem(1)); % Expected number of full elements
nup = done-ntop; % Number of elements that still need to be filled
xnew = zeros(length(x),1); 

if (nup<=0) % Check if the number of filled elements is greater or equal to the limit
    [~,su] = sort(x(xtop),'descend');
    xnew(xtop(su(1:done))) = 1.0; % The largest "done" densities are rounded to 1
    return;
else
    xnew(xtop) = 1.0; % Only the densities that are greater or equal to "ut" are rounded to 1
end

xvar = find((x>prj.lt)&(x<prj.ut)); % Indexes of the intermediate densities
nvar = length(xvar);
if (nvar==0) % Check if all of the densities are sufficiently close to 0 or 1
    return
end
dvol = min([nvar,nup,round(volfrac*length(x))-ntop]); % Number of element densities that need to be rounded to 1

if (prj.nearlim) % Checking if the densities are to be rounded to the nearest limit (0 or 1)
    
    % Trying to round to one the largest volfrac*length(x) variables and to zero the remaining variables
    [~,ivol] = sort(x(xvar),'descend');
    xnew(xvar(ivol(1:dvol))) = 1.0;
    
    % Verifying if we have a descent direction
    d = xnew - x;
    cosang = d'*g/(gnorm*norm(d));
else
    cosang = -1.0; % Dummy value for the cosine of the angle between d and -grad(Lag)
end

if (cosang < coslim) % Checking if we still don't have a descent direction

    % Trying to move into the direction defined by the projection of -grad(Lag) onto the box

%     iup = find(g(xvar)>gt); % Indexes of the intermediate densities that can be rounded up to 1
%     iup2 = find((abs(g(xvar))<=gt)&(x(xvar)>=prj.vumin)); % Indexes of the densities with g = 0 that can be rounded up to 1
    iup = find(g(xvar)>0.0); % Indexes of the intermediate densities that can be rounded up to 1
    iup2 = find((g(xvar)==0.0)&(x(xvar)>=prj.vumin)); % Indexes of the densities with g = 0 that can be rounded up to 1
    g(xvar(iup2)) = 1e16; % The gradient of xvar(iup2) is set to 1e16 so these variables go to the front of the queue
    iup = [iup;iup2]; % Complete list of densities that can be rounded up to 1
    niup = length(iup);
    idown = find(g(xvar)<-gt);
    nidown = length(idown);

    alphau = (1-x(xvar(iup)))./g(xvar(iup)); % Step lengths for the densities that can be rounded up to 1
    alphad = -x(xvar(idown))./g(xvar(idown)); % Step lengths for the densities that can be rounded down to 0
    [alphau,su] = sort(alphau);
    [alphad,sd] = sort(alphad);
    xup = x(xvar(iup(su))); % Ordered list of the densities that can be rounded up to 1
    xdown = x(xvar(idown(sd))); % Ordered list of the densities that can be rounded down to 0
    kup = find(xup<prj.vumin,1); % First density in xup that violates the condition xup>=vumin
    kdown = find(xdown>prj.vlmax,1); % First density in xdown that violates the condition xdown<=vlmax

    if (isempty(kup))
        kup = niup; % All of the densities can be rounded up to 1
    else
        kup = kup - 1; % kup is the index of the last density that can be rounded up to 1
    end
    if isempty(kdown)
        kdown = nidown; % All of the densities can be rounded down to 0
    else
        kdown = kdown - 1; % kdown is the index of the last density that can be rounded down to 0
    end

    if (niup==0) 
        if (nidown==0)
            return % There are no densities that can be rounded to 0 or 1
        else
            alpha = alphad(1:kdown); % All of the densities will be rounded down to 0
        end
    else
        imax = min(dvol,kup); % Actual number of densities that can be rounded up to 1
        if (nidown==0)
            alpha = alphau(1:imax); % All of the densities will be rounded up to 1
        elseif (imax==0)
            alpha = alphad(1:kdown); % All of the densities will be rounded down to 0
        else
            imax2 = find(alphad>alphau(imax),1); 
            if isempty(imax2)
                imax2 = kdown; % There are kdown densities that can be rounded down to 0
            else
                imax2 = min(kdown,imax2-1); % imax2 is the index of the last density that will be rounded down to 0
            end
            alpha = sort([alphau(1:imax);alphad(1:imax2)]); % Ordered list of step lengths
        end
    end

    pos = length(alpha);
    while ((pos>0)&&(cosang<coslim))

        % Trying to take the largest step that satifies the angle condition

        xnew(xvar) = max(0.0,min(1.0,x(xvar)+alpha(pos)*g(xvar)));
        d = xnew - x;
        cosang = d'*g/(gnorm*norm(d));
        pos = pos-1;
    end

    if (cosang<coslim) % Checking if all of our attempts have failed
        xnew(xvar) = x(xvar); % Only the densities thar are greater than ut or smaller than lt are rounded.
    end

end

end