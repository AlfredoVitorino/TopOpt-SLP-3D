function Display3D(str,rho,prjop)
% Display3D creates a 3D view of the optimal structure using the element densities vector. %
% INPUT: str - structure with the problem data.
%        rho - element densities vector.
%        prjop - indicates if the densities were rounded to 0 or 1 at the end of the SLP algorithm.
% ---------- %

ny = str.ny; 
nxny = str.nx*str.ny;

facesXY = zeros(2*str.nelem,4); 
facesXYcolor = zeros(2*str.nelem,3);
j = 1;
for i = 1:str.nelem 
    facesXY(i,1) = j;
    if(mod(j,ny) == (ny-1))
        j = j+2;
    else
        j = j+1;
    end   
    if(mod(j,nxny) == (nxny-str.nely))
        j = j+ny;
    end
end 

if (prjop)
    color = 0.25+0.75*(1-rho); 
else
    color = 0.1+0.9*(1-rho);
end
facesXYcolor(1:str.nelem,:) = [color, color, color]; 

facesYZ = facesXY; 
facesZX = facesXY; 
facesYZcolor = facesXYcolor; 
facesZXcolor = facesXYcolor;  

facesXY(str.nelem+1:end,1) = facesXY(1:str.nelem,1)+nxny;
facesXYcolor(str.nelem+1:end,:) = facesXYcolor(1:str.nelem,:);
facesXY(:,2:end) = [facesXY(:,1)+ny, facesXY(:,1)+ny+1, facesXY(:,1)+1];
rhoVis = facesXYcolor(:,1) <= 0.5; 
facesXY = facesXY(rhoVis,:); 
facesXYcolor = facesXYcolor(rhoVis,:);

facesYZ(str.nelem+1:end,1) = facesYZ(1:str.nelem,1)+ny;
facesYZcolor(str.nelem+1:end,:) = facesYZcolor(1:str.nelem,:);
facesYZ(:,2:end) = [facesYZ(:,1)+1, facesYZ(:,1)+nxny+1, facesYZ(:,1)+nxny];
rhoVis = facesYZcolor(:,1) <= 0.5; 
facesYZ = facesYZ(rhoVis,:); 
facesYZcolor = facesYZcolor(rhoVis,:);

facesZX(str.nelem+1:end,1) = facesZX(1:str.nelem,1)+1; 
facesZXcolor(str.nelem+1:end,:) = facesZXcolor(1:str.nelem,:);
facesZX(:,2:end) = [facesZX(:,1)+ny, facesZX(:,1)+nxny+ny, facesZX(:,1)+nxny];
rhoVis = facesZXcolor(:,1) <= 0.5; 
facesZX = facesZX(rhoVis,:); 
facesZXcolor = facesZXcolor(rhoVis,:);

figure

[X,Y,Z] = meshgrid(0:str.nelx,0:str.nely,0:str.nelz);
patch('Faces', [facesXY; facesYZ; facesZX], 'Vertices', [X(:), Z(:), Y(:)], 'FaceColor', 'flat', 'FaceVertexCData', [facesXYcolor; facesYZcolor; facesZXcolor]);

set(gca,'Ydir','reverse');
axis tight; axis equal; axis off; box on; rotate3d on; 
view(30,20);

end