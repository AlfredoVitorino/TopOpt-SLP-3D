function Display3D(str,rho,prjop)
% Display3D creates a 3D view of the optimal structure using the element densities vector %
% INPUT: str - structure with the problem data.
%        rho - element densities vector.
%        prjop - indicates if the densities were rounded to 0 or 1 at the end of the SLP algorithm.
% ---------- %

elemvis = zeros(str.nelem, 1); % indicates the visible elements (the ones solid and not surrounded by solid elements) 
ex = 1; ey = 1; ez = 1; % layers of the element in the x, y and z directions on the mesh
for e = 1:str.nelem   
    if(rho(e) >= 0.5)
        ex2 = [ex, ex-1, ex, ex, ex+1, ex]; 
        ey2 = [ey, ey, ey-1, ey+1, ey, ey]; 
        ez2 = [ez-1, ez, ez, ez, ez, ez+1];  
        
        count = 0; 
        for j = 1:6
            i2 = ex2(j); 
            j2 = ey2(j); 
            k2 = ez2(j); 
            if(i2 >= 1 && i2 <= str.nelx && j2 >= 1 && j2 <= str.nely && k2 >= 1 && k2 <= str.nelz) 
                e2 = (k2-1)*str.nelx*str.nely + (i2-1)*str.nely + j2;
                if(rho(e2) >= 0.5)
                    count = count + 1; 
                end
            end
        end     
        
        if(count < 6)
            elemvis(e) = 1; 
        end
    end
    
    if(mod(e,str.nelx*str.nely) == 0)
        ex = 1;
        ey = 1;
        ez = ez + 1; 
    elseif(mod(e,str.nely) == 0)
        ex = ex + 1;
        ey = 1; 
    else
        ey = ey + 1; 
    end
end
rho(elemvis == 0) = 0; % consider that all not visible elements have zero density 

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
    %color = 0.25+0.75*(1-rho); 
    color = 0.45+0.55*(1-rho); 
else
    color = 0.1+0.9*(1-rho);
end
facesXYcolor(1:str.nelem,:) = [color, color, color]; 

facesYZ = facesXY; 
facesZX = facesXY; 
facesYZcolor = facesXYcolor; 
facesZXcolor = facesXYcolor;  

rho2 = [rho;rho]; 

facesXY(str.nelem+1:end,1) = facesXY(1:str.nelem,1)+nxny;
facesXYcolor(str.nelem+1:end,:) = facesXYcolor(1:str.nelem,:);
facesXY(:,2:end) = [facesXY(:,1)+ny, facesXY(:,1)+ny+1, facesXY(:,1)+1];
%rhoVis = facesXYcolor(:,1) <= 0.5; 
rhoVis = rho2 >= 0.5; 
facesXY = facesXY(rhoVis,:); 
facesXYcolor = facesXYcolor(rhoVis,:);

facesYZ(str.nelem+1:end,1) = facesYZ(1:str.nelem,1)+ny;
facesYZcolor(str.nelem+1:end,:) = facesYZcolor(1:str.nelem,:);
facesYZ(:,2:end) = [facesYZ(:,1)+1, facesYZ(:,1)+nxny+1, facesYZ(:,1)+nxny];
%rhoVis = facesYZcolor(:,1) <= 0.5; 
rhoVis = rho2 >= 0.5; 
facesYZ = facesYZ(rhoVis,:); 
facesYZcolor = facesYZcolor(rhoVis,:);

facesZX(str.nelem+1:end,1) = facesZX(1:str.nelem,1)+1; 
facesZXcolor(str.nelem+1:end,:) = facesZXcolor(1:str.nelem,:);
facesZX(:,2:end) = [facesZX(:,1)+ny, facesZX(:,1)+nxny+ny, facesZX(:,1)+nxny];
%rhoVis = facesZXcolor(:,1) <= 0.5; 
rhoVis = rho2 >= 0.5; 
facesZX = facesZX(rhoVis,:); 
facesZXcolor = facesZXcolor(rhoVis,:);

figure

[X,Y,Z] = meshgrid(0:str.nelx,0:str.nely,0:str.nelz);
%patch('Faces', [facesXY; facesYZ; facesZX], 'Vertices', [X(:), Z(:), Y(:)], 'FaceColor', 'flat', 'FaceVertexCData', [facesXYcolor; facesYZcolor; facesZXcolor]);
patch('Faces', [facesXY; facesYZ; facesZX], 'Vertices', [X(:), Z(:), Y(:)], 'FaceColor', 'flat', 'LineWidth', 0.001, 'EdgeLighting', 'flat', 'FaceVertexCData', [facesXYcolor; facesYZcolor; facesZXcolor]);

set(gca,'Ydir','reverse');
axis tight; axis equal; axis off; box on; rotate3d on; 
view(30,20);
 
light('Position', [-1,0.8,0]); 
light('Position', [-1,0,0.8]);

end