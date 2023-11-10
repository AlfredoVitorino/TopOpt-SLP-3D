function [elem2fix0,elem2fix1,dofs2fix,dsgn2fix,dens2fix] = Elem2Fix(str,strDsgn,strDens,x,xfil,mr,Dofs,genK,elem,Fgrad)
% Elem2Fix choose the elements that will have fixed variables when solving the problem again increasing the element degree. %
% INPUT: str - structure with the problem data.
%        strDens - structure with density mesh data for multiresolution. 
%        strDsgn - structure with design variable mesh data for multiresolution.
%        x - design variables vector.
%        xfil - element densities vector.
%        mr - structure that contains parameters used by the multiresolution method.
%        Dofs - matrix with the indexes of the degrees of freedom for each element.
%        genK - explicit generation of the global stiffness matrix (true = yes, false = no).
%        elem - structure with element characteristics.
%        Fgrad - gradient of the objective function.
% OUTPUT: elem2fix0 - void displacement elements that will be fixed.
%         elem2fix1 - solid displacement elements that will be fixed.
%         dofs2fix - degrees of freedom that will be eliminated from the equilibrium linear system. 
%         dsgn2fix - design variables that will be fixed. 
%         dens2fix - density elements that will be fixed. 
% ---------- % 

elem2fix = zeros(str.nelem,1); % displacement elements that will be fixed
elem2fix0 = zeros(str.nelem,1); % void displacement elements that will be fixed
elem2fix1 = zeros(str.nelem,1); % solid displacement elements that will be fixed
dofs2fix = zeros(3*str.nnodes,1); % degrees of freedom that will be fixed
dsgn2fix = zeros(strDsgn.nelem,1); % design variables that will be fixed
dens2fix = zeros(strDens.nelem,1); % density elements that will be fixed

DsgnInd = zeros(str.nelem, mr.d^3); % indexes of the design variables inside each displacement element
DensInd = zeros(str.nelem, mr.n^3); % indexes of the densities inside each displacement element

nelyd = strDsgn.nely; 
nelxyd = strDsgn.nely*strDsgn.nelx;
nelyn = strDens.nely; 
nelxyn = strDens.nely*strDens.nelx;

% Find the displacement elements that are full of approximately zero (or one) densities % 
fullelem = zeros(str.nelem,1); % indicates if the displacement element is void (=1) or solid (=2)
elemLayer = zeros(str.nelem,3); % layers of the elements on the displacement mesh
ex = 1; ey = 1; ez = 1;
for i = 1:str.nelem 
    % Find the indexes of the design variables inside the displacement element i
    vbd = zeros(mr.d^2, 1); 
    indd = 1; 
    for n1 = 0:(mr.d-1) 
        for n2 = 0:(mr.d-1)
            vbd(indd) = strDsgn.fdv(i) + n1*nelxyd + n2*nelyd; 
            indd = indd + 1; 
        end
    end         
    vd = zeros(mr.d^2, mr.d);
    for n3 = 0:(mr.d-1)
        vd(:,n3+1) = vbd+n3; 
    end
    vd = vd(:);
    DsgnInd(i,:) = vd; 
    
    % Find the indexes of the density elements inside the displacement element i   
    vbn = zeros(mr.n^2, 1); 
    indn = 1; 
    for n1 = 0:(mr.n-1) 
        for n2 = 0:(mr.n-1)
            vbn(indn) = strDens.fde(i) + n1*nelxyn + n2*nelyn; 
            indn = indn + 1; 
        end
    end         
    vn = zeros(mr.n^2, mr.n);
    for n3 = 0:(mr.n-1)
        vn(:,n3+1) = vbn+n3; 
    end
    vn = vn(:); 
    DensInd(i,:) = vn; 
     
    % Check if the element is void or solid 
    if(elem.fixdsgn)
        maxdens = max(x(vd));
        mindens = min(x(vd));
    else
        maxdens = max(xfil(vn));
        mindens = min(xfil(vn));
    end
    if(maxdens <= elem.fixtl && min(Fgrad(vd)) >= -elem.tolgdsgn) % the element is void and the design variables shall not increase according to the gradient
        fullelem(i) = 1; 
    elseif(mindens >= elem.fixtu && max(Fgrad(vd)) <= elem.tolgdsgn) % the element is solid and the design variables shall not decrease according to the gradient
        fullelem(i) = 2; 
    end
    
    % Update the layers of the displacement element 
    elemLayer(i,:) = [ex ey ez];   
    if(mod(i,str.nelx*str.nely) == 0)
        ex = 1;
        ey = 1;
        ez = ez + 1; 
    elseif(mod(i,str.nely) == 0)
        ex = ex + 1;
        ey = 1; 
    else
        ey = ey + 1; 
    end
end
% ---------- % 

% Find the void (solid) elements that are surrounded by void (solid) elements %
elm = 1:str.nelem;
if(~elem.fixnb)
    for e1 = elm(fullelem ~= 0)
        i1 = elemLayer(e1,1);
        j1 = elemLayer(e1,2);
        k1 = elemLayer(e1,3);

        cont = true;
        k2 = max(k1-1, 1);
        kend = min(k1+1, str.nelz);
        while ((k2<=kend)&&cont)
            i2 = max(i1-1, 1);
            iend = min(i1+1, str.nelx);
            while ((i2<=iend)&&cont)
                j2 = max(j1-1, 1);
                jend = min(j1+1, str.nely);
                while ((j2<=jend)&&cont)
                    e2 = (k2-1)*str.nelx*str.nely + (i2-1)*str.nely + j2;
                    if(fullelem(e2) ~= fullelem(e1))
                        if (~elem.fixdsgn)
                            cont = false; 
                        elseif (fullelem(e1) == 1)
                            de2 = DsgnInd(e2,:);
                            if ((fullelem(e2) == 2)||(min(Fgrad(de2))<-elem.tolgdsgn)||(max(x(de2))>elem.ltdsgn))
                                cont = false; 
                            end
                        else % In this case, fullelem(e1) == 2
                            de2 = DsgnInd(e2,:);
                            if ((fullelem(e2) == 1)||(max(Fgrad(de2))>elem.tolgdsgn)||(min(x(de2))<elem.utdsgn))
                                cont = false; 
                            end
                        end
                    end
                    j2 = j2+1;
                end
                i2 = i2+1;
            end
            k2 = k2+1;
        end
        if(cont)
            elem2fix(e1) = 1;
        end
    end
else % Choose all void and solid elements to be fixed
    elem2fix(fullelem ~= 0) = 1;
end
% ---------- % 

% Set the design variables and densities that will be fixed % 
for i = elm(elem2fix == 1)  
    dsgnind = DsgnInd(i,:);
    densind = DensInd(i,:);     
    dsgn2fix(dsgnind) = 1; 
    dens2fix(densind) = 1;
    if(fullelem(i) == 1)
        elem2fix0(i) = 1; 
    elseif(fullelem(i) == 2)
        elem2fix1(i) = 1;
    end
end
% ---------- %

% Set the degrees of freedom that will be eliminated from the linear system %
if(elem.fixdofs ~= 0 && genK)
    if(elem.fixdofs == 1) 
        dofs2fix = ones(3*str.nnodes,1);
    end
    for i = 1:str.nelem
        if(elem.fixdofs == 1 && (elem2fix(i) == 0 || fullelem(i) == 2)) % Do not fix the displacements of solid elements or nodes shared by neighbor elements 
            dofs2fix(Dofs(i,:)) = 0;
        elseif(elem.fixdofs == 2 && (elem2fix(i) == 1 && fullelem(i) == 1)) % Fix all displacements of the void elements        
            dofs2fix(Dofs(i,:)) = 1;
        end
    end
end
% ---------- % 

% % Use this code if you want to see a figure with the displacement elements that are fully solid or void % 
% Display3D(strDens,xfil,false);
% hold on;
% face = [1,2,3,4; 5,6,7,8; 1,5,8,4; 2,6,7,3; 1,2,6,5; 4,3,7,8];
% ex = 0; ey = 0; ez = 0; j = 1;
% for i = 1:str.nelem 
%     if(elem2fix(i) == 1)   
%         if(fullelem(i) == 1)
%             vert = (mr.n)*[ex,ey,ez; ex+str.el,ey,ez; ex+str.el,ey+str.eh,ez; ex,ey+str.eh,ez; ex,ey,ez+str.ew; ex+str.el,ey,ez+str.ew; ex+str.el,ey+str.eh,ez+str.ew; ex,ey+str.eh,ez+str.ew]; 
%             vert(:,[2,3]) = vert(:,[3,2]);
%             patch('Faces', face, 'Vertices', vert, 'FaceColor', [0, 0, 1], 'EdgeColor', [0, 0, 1], 'FaceAlpha', 0.5, 'LineWidth',2); 
%             hold on;
%             plot3((mr.n)*(ex+str.el/2),(mr.n)*(ez+str.ew),(mr.n)*(ey+str.eh/2),'bo'); 
%             hold on;
%         elseif(fullelem(i) == 2)
%             vert = (mr.n)*[ex,ey,ez; ex+str.el,ey,ez; ex+str.el,ey+str.eh,ez; ex,ey+str.eh,ez; ex,ey,ez+str.ew; ex+str.el,ey,ez+str.ew; ex+str.el,ey+str.eh,ez+str.ew; ex,ey+str.eh,ez+str.ew]; 
%             vert(:,[2,3]) = vert(:,[3,2]);
%             patch('Faces', face, 'Vertices', vert, 'FaceColor', [1, 0, 0], 'EdgeColor', [1, 0, 0], 'FaceAlpha', 0.5, 'LineWidth',2); 
%             hold on;
%             plot3((mr.n)*(ex+str.el/2),(mr.n)*(ez+str.ew),(mr.n)*(ey+str.eh/2),'rx'); 
%             hold on;
%         end
%     end
%     if(mod(j,str.nely+1) == str.nely)
%         ex = ex+str.el; ey = 0; j = j+2;
%     else
%         ey = ey+str.eh; j = j+1;
%     end   
%     if(mod(j,(str.nelx+1)*(str.nely+1)) == ((str.nelx+1)*(str.nely+1)-str.nely))
%         ex = 0; ey = 0; ez = ez+str.ew; j = j+str.nely+1;
%     end
% end 
% set(gca,'Ydir','reverse');
% axis tight; axis equal; axis off; box on; rotate3d on; 
% view(30,20);
% 
% % Use this code if you want to see a figure with the density elements that were choosen to be fixed %
% Display3D(strDens,xfil,false);
% hold on;
% face = [1,2,3,4; 5,6,7,8; 1,5,8,4; 2,6,7,3; 1,2,6,5; 4,3,7,8];
% ex = 0; ey = 0; ez = 0; j = 1;
% for i = 1:str.nelem 
%     if(elem2fix(i) == 1)
%         vert = (mr.n)*[ex,ey,ez; ex+str.el,ey,ez; ex+str.el,ey+str.eh,ez; ex,ey+str.eh,ez; ex,ey,ez+str.ew; ex+str.el,ey,ez+str.ew; ex+str.el,ey+str.eh,ez+str.ew; ex,ey+str.eh,ez+str.ew]; 
%         vert(:,[2,3]) = vert(:,[3,2]);
%         patch('Faces', face, 'Vertices', vert, 'FaceColor', [0.9 0.3 0.1], 'FaceAlpha', 0.5); 
%         hold on;
%         xc = (mr.n)*[ex+(1/4)*str.el, ex+(3/4)*str.el, ex+(1/4)*str.el, ex+(3/4)*str.el]; 
%         yc = (mr.n)*[ez+str.ew, ez+str.ew, ez+str.ew, ez+str.ew]; 
%         zc = (mr.n)*[ey+(1/4)*str.eh, ey+(1/4)*str.eh, ey+(3/4)*str.eh, ey+(3/4)*str.eh]; 
%         plot3(xc,yc,zc,'*','Color',[0.9 0.3 0.1]); 
%         hold on;
%     end
%     if(mod(j,str.nely+1) == str.nely)
%         ex = ex+str.el; ey = 0; j = j+2;
%     else
%         ey = ey+str.eh; j = j+1;
%     end   
%     if(mod(j,(str.nelx+1)*(str.nely+1)) == ((str.nelx+1)*(str.nely+1)-str.nely))
%         ex = 0; ey = 0; ez = ez+str.ew; j = j+str.nely+1;
%     end
% end 
% set(gca,'Ydir','reverse');
% axis tight; axis equal; axis off; box on; rotate3d on; 
% view(30,20);
% 
% % Use this code if you want to see a figure with the nodes whose degrees of freedom were choosen to be fixed (works only for Lagrange elements with degree 2) %
% Display3D(strDens,xfil,false);
% hold on;
% face = [1,2,3,4; 5,6,7,8; 1,5,8,4; 2,6,7,3; 1,2,6,5; 4,3,7,8];
% ex = 0; ey = 0; ez = 0; j = 1;
% xc = []; yc = []; zc = [];
% xnc = []; ync = []; znc = []; 
% for i = 1:str.nelem 
%     vert = (mr.n)*[ex,ey,ez; ex+str.el,ey,ez; ex+str.el,ey+str.eh,ez; ex,ey+str.eh,ez; ex,ey,ez+str.ew; ex+str.el,ey,ez+str.ew; ex+str.el,ey+str.eh,ez+str.ew; ex,ey+str.eh,ez+str.ew]; 
%     vert(:,[2,3]) = vert(:,[3,2]);
%     patch('Faces', face, 'Vertices', vert, 'FaceColor', 'none');
%     if(fullelem(i) == 1 && elem2fix(i) == 1)
%         xc = [xc; (mr.n)*repmat([ex; ex+(1/2)*str.el; ex+str.el],9,1)]; 
%         yc = [yc; (mr.n)*[repmat(ez,9,1); repmat(ez+(1/2)*str.ew,9,1); repmat(ez+str.ew,9,1)]]; 
%         zc = [zc; (mr.n)*repmat([ey; ey; ey; ey+(1/2)*str.eh; ey+(1/2)*str.eh; ey+(1/2)*str.eh; ey+str.eh; ey+str.eh; ey+str.eh],3,1)]; 
%     elseif(elem.fixdofs == 1)
%         xnc = [xnc; (mr.n)*repmat([ex; ex+(1/2)*str.el; ex+str.el],9,1)]; 
%         ync = [ync; (mr.n)*[repmat(ez,9,1); repmat(ez+(1/2)*str.ew,9,1); repmat(ez+str.ew,9,1)]]; 
%         znc = [znc; (mr.n)*repmat([ey; ey; ey; ey+(1/2)*str.eh; ey+(1/2)*str.eh; ey+(1/2)*str.eh; ey+str.eh; ey+str.eh; ey+str.eh],3,1)];
%     end
%     if(mod(j,str.nely+1) == str.nely)
%         ex = ex+str.el; ey = 0; j = j+2;
%     else
%         ey = ey+str.eh; j = j+1;
%     end   
%     if(mod(j,(str.nelx+1)*(str.nely+1)) == ((str.nelx+1)*(str.nely+1)-str.nely))
%         ex = 0; ey = 0; ez = ez+str.ew; j = j+str.nely+1;
%     end
% end 
% c = [xc,yc,zc];
% if(elem.fixdofs == 1)
%     nc = [xnc,ync,znc];
%     ind = ismember(c,nc,'rows');
%     c(ind==1,:) = [];
% end
% xc = c(:,1); yc = c(:,2); zc = c(:,3); 
% plot3(xc,yc,zc,'*','Color',[0 0.7 0]);  
% set(gca,'Ydir','reverse');
% axis tight; axis equal; axis off; box on; rotate3d on; 
% view(30,20);

elm = 1:str.nelem;
elem2fix0 = elm(elem2fix0 == 1); 
elem2fix1 = elm(elem2fix1 == 1); 
dofs = 1:(3*str.nnodes);
dofs2fix = dofs(dofs2fix == 1);
dsgn = 1:strDsgn.nelem;
dsgn2fix = dsgn(dsgn2fix == 1); 
dens = 1:strDens.nelem;
dens2fix = dens(dens2fix == 1);

disp(['Fixed dofs: ',num2str(length(dofs2fix)), '. Fixed design variables: ',num2str(length(dsgn2fix)), '. Fixed densities: ',num2str(length(dens2fix))]); 

end