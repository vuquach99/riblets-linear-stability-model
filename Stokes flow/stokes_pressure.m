%% Programme to solve for Stokes velocity profile in riblets
% Solve A*u = P
% S = riblets grid points
% P = pressure gradient
% A = Laplacian coefficients
% Sc = Sd - S = riblets surface grid points

clear
close all
tic

%% Input data
% Sides of the non-periodic-in-z domain
ly = 0.010186218918508; % always <= height for pressure-driven flow
% blade: 0.368504583707110
% triangle9: 0.351838677288411
% triangle6: 0.690839905195580
% triangle3: 1.666214643830251
% trapezium: 0.292484372828901

lz = 1; % always

% Define grid sizes
n = 100;
nz = 2*n+1; % number of points in non-periodic z (spanwise)
ny = round(ly*(nz-1)); % number of points y (wall-normal) domain is enclosed
dy = ly/(ny-1); % length of sub-intervals in y-axis
dz = lz/(nz-1); % length of sub-intervals in z-axis
z = linspace(0,lz,nz);
y = linspace(0,ly,ny);

dpdx = -1; % Pressure gradient
Sx = 0; % Shear at the top, could set to 1 (if normalised, always 1)

geometry = 3; 
% 0 = parabola (k/s = 1)
% 1 = triangle (k/s = 1.866 for 30deg, 0.5 for 90deg, 0.866 for 60deg)
% 2 = semi-circle (k/s = 0.5)
% 3 = trapezium (k/s = 0.5; tip half-angle = 15deg)
% 4 = blade (k/s = 0.5; t/s = 0.2)
angle = 60; % 30/60/90 degrees, for triangles only

savefile = 0;

%% Build S matrix - Grid Points of riblets
% S = 1 for points within and on boundary and = 0 elsewhere 
S = zeros(ny,nz);
if geometry == 0 % parabola
    parabola = (2*(z-lz/2)).^2;
    shape = 'parabola';
    height = 1;
    for k = 1:nz
        for j = 1:ny
            if y(j) < parabola(k)
                S(j,k) = 1;
            else
                S(j,k) = 0;
            end
        end
    end
end

if geometry == 1 % triangle
    if angle == 90        
        shape = 'triangle9';
        height = 0.5;
        triangle = z-lz/2;
        triangle2 = -z+lz/2;
    end
    if angle == 60
        shape = 'triangle6';
        height = sqrt(3)/2;
        triangle = sqrt(3)*z-sqrt(3)*lz/2;
        triangle2 = -sqrt(3)*z+sqrt(3)*lz/2;
    end
    if angle == 30
        shape = 'triangle3';
        height = 1 + sqrt(3)/2;
        triangle = (2+sqrt(3))*z-1-sqrt(3)*lz/2;
        triangle2 = -(2+sqrt(3))*z+1+sqrt(3)*lz/2;
    end
        for k = 1:nz
            for j = 1:ny
                if y(j) < triangle(k) || y(j) < triangle2(k)
                    S(j,k) = 1;
                else
                    S(j,k) = 0;
                end
            end
        end
end

if geometry == 2 % semi-circle
    shape = 'circle';
    height = 0.5;
    circle = (-sqrt((lz/2)^2-(z-lz/2).^2))+lz/2;
    for k = 1:nz
        for j = 1:ny
            if y(j) < circle(k)
                S(j,k) = 1;
            else
                S(j,k) = 0;
            end
        end
    end
end

if geometry == 3 % trapezium
    shape='trapezium';
    height = 0.5;
    trapezium = (4+2*sqrt(3))*z-3.5-2*sqrt(3);
    trapezium2 = -(4+2*sqrt(3))*z+0.5;
    for k=1:nz
        for j=1:ny
            if y(j) < trapezium(k) || y(j) < trapezium2(k)
                S(j, k) = 1;
            else
                S(j, k) = 0;
            end
        end
    end
end

if geometry == 4 % blade
    shape = 'blade';
    height = 0.5;
    for j = 1:height*nz
        for k = 1:nz
            if z(k) < 0.1*lz || z(k) > 0.9*lz
                S(j,k) = 1;
            else
                S(j,k) = 0;
            end
        end
    end
end

figure
spy(S)
title('S matrix')

%% Build Sd matrix
% Sd = 1 for points within and on boundary and = 0 elsewhere 
% riblets now bigger by rd = max(dy,dz)
rd = max(dy,dz);
Sd = zeros(ny,nz);
if geometry == 0 % parabola
    parabola_d = (2*(z-lz/2)).^2+rd;
    for k = 1:nz
        for j = 1:ny
            if y(j) < parabola_d(k)
                Sd(j,k) = 1;
            else
                Sd(j,k) = 0;
            end
        end
    end
end

if geometry == 1 % triangle
    if angle == 90
        triangle = z-lz/2+rd;
        triangle2 = -z+lz/2+rd;
    end
    if angle == 60
        triangle = sqrt(3)*z-sqrt(3)*lz/2+rd;
        triangle2 = -sqrt(3)*z+sqrt(3)*lz/2+rd;
    end
    if angle == 30
        triangle = (2+sqrt(3))*z-1-sqrt(3)*lz/2+rd;
        triangle2 = -(2+sqrt(3))*z+1+sqrt(3)*lz/2+rd;
    end
    for k = 1:nz
        for j = 1:ny
            if y(j)<triangle(k) || y(j)<triangle2(k)
                Sd(j,k) = 1;
            else
                Sd(j,k) = 0;
            end
        end
    end
end

if geometry == 2 % semi-circle
    circle = (-sqrt((lz/2)^2-(z-lz/2).^2))+lz/2+rd;
    for k = 1:nz
        for j = 1:ny
            if y(j)<circle(k)
                Sd(j,k) = 1;
            else
                Sd(j,k) = 0;
            end
        end
    end
end

if geometry == 3 % trapezium
    trapezium = (4+2*sqrt(3))*z-3.5-2*sqrt(3)+rd;
    trapezium2 = -(4+2*sqrt(3))*z+0.5+rd;
    for k=1:nz
        for j=1:ny
            if y(j)<trapezium(k) || y(j)<trapezium2(k)
                Sd(j, k) = 1;
            else
                Sd(j, k) = 0;
            end
        end
    end
    Sd(1,:) = ones(1,nz);
end

if geometry == 4 % blade
    for j = 1:height*nz %+1 if domain height > 0.5*lz
        for k = 1:nz
            if z(k) < 0.1*lz+dz-0.00001 || z(k) > 0.9*lz-dz+0.00001
                Sd(j,k) = 1;
            else
                Sd(j,k) = 0;
            end
        end
    end
    Sd(1,:) = ones(1,nz);
end

figure
spy(Sd)
title('Sd matrix')

%% Build Sc matrix - curve of the boundary
Sc = Sd-S;
figure
spy(Sc)
title('Sc matrix')

%% Build P matrix - pressure gradient - RHS of Stokes
P = ones(ny*nz,1);
P = sparse(P);
P = P*dpdx;

%% Build A matrix 
beta = -2/(dy^2)-2/(dz^2); % coefficient of point p
alpha = 1/(dy^2); % coefficient of surrounding points
r1 = beta.*ones(ny,1);
r2 = alpha.*ones(ny-1,1);
r3 = alpha.*ones(1,1);
B = diag(r1,0) + diag(r2,1) + diag(r2,-1) + diag(r3,-ny+1) + diag(r3,ny-1);
B = sparse(B);
C = kron(eye(nz),B);
C = sparse(C);
gamma = 1/(dz^2);
r4 = gamma.*ones(ny*nz-ny,1);
r5 = gamma.*ones(ny,1);
D = diag(r4,-ny)+diag(r4,+ny)+diag(r5,-ny*nz+ny)+diag(r5,ny*nz-ny);
D = sparse(D);
A = C+D;
A = sparse(A);

%% Shear at the top condition
for k = 1:nz
    j = ny; % top of domain
    index = (k-1)*ny+j;
    A(index,:) = 0;
    A(index,index) = 1/dy^2; % should be dy, but P(index) changed as well
    if index~=1
        A(index,index-1) = -1/dy^2;
    end   
    P(index) = Sx/dy;
end

%% Solid condition - modify the A matrix
% Modify A and P matrix for points within and on the boundary of riblets
% Velocity inside riblets = 0
for j = 1:ny
    for k = 1:nz
        if S(j,k) == 1
            index = (k-1)*ny+j;
            A(index,:) = 0;
            A(index,index) = 1/dy^2;
            P(index) = 0;
        end
    end
end

%% Modify the A matrix - curve of the boundary
% Modification of the matrix A near the boundary of the solid
for j=1:ny
    for k=1:nz
        if Sc(j,k) == 1
            y0 = (j-ny)*dy;
            z0 = k*dz;          
            my = abs(floor(log10(dy)));
            mz = abs(floor(log10(dz)));
            
        %ys: (j-+2,k);(j-+1,k);(j,k-+2);(j,k-+1) 
            ys=[y0-dy, y0, y0+dy];
            zs=[z0-dz, z0, z0+dz];  
            
            if geometry == 0 % parabola
                yp0 = (2*(z0-lz/2)).^2;
                zpPlus = sqrt(y0)/2+lz/2;
                zpMinus = -sqrt(y0)/2+lz/2;
            end
            
            if geometry == 1 % triangle
                if angle == 90
                    if z0 <= lz/2
                        yp0 = -z0+lz/2;
                    else
                        yp0 = z0-lz/2;
                    end
                    zpPlus = y0+lz/2;
                    zpMinus = -y0+lz/2;
                end
                if angle == 60
                    if z0 <= lz/2
                        yp0 = -sqrt(3)*z0+sqrt(3)*lz/2;
                    else
                        yp0 = sqrt(3)*z0-sqrt(3)*lz/2;
                    end
                    zpPlus = y0/sqrt(3)+lz/2;
                    zpMinus = -y0/sqrt(3)+lz/2;
                end
                if angle == 30
                    if z0 <= lz/2
                        yp0 = (2+sqrt(3))*z0-1-sqrt(3)*lz/2;
                    else
                        yp0 = -(2+sqrt(3))*z0+1+sqrt(3)*lz/2;
                    end
                    zpPlus = y0/(2+sqrt(3))+(1+sqrt(3)*lz/2)/(2+sqrt(3));
                    zpMinus = -y0/(2+sqrt(3))+(1+sqrt(3)*lz/2)/(2+sqrt(3));
                end
            end
            
            if geometry == 2 % semi-circle
                yp0=-(-sqrt((lz/2)^2-(z0-lz/2).^2))+ly;
                zpPlus=sqrt((lz/2)^2-(y0-ly).^2)+lz/2;
                zpMinus=-zpPlus;
            end
            
            if geometry == 3 % trapezium
                if z0 <= 0.5-sqrt(3)/4
                    yp0 = -(4+2*sqrt(3))*z+0.5;
                elseif z0 >= 0.5+sqrt(3)/4
                    yp0 = (4+2*sqrt(3))*z-3.5-2*sqrt(3);
                else
                    yp0 = 0;
                end
                zpPlus = (3.5+2*sqrt(3)+y0)/(4+2*sqrt(3));
                zpMinus = (0.5-y0)/(4+2*sqrt(3));
            end
            
            if geometry == 4 % blade
                if z0 < 0.1*lz || z0 > 0.9*lz
                    yp0 = height*lz;
                else
                    yp0 = 0;
                end
                if y0 <= height*lz
                    zpPlus = 0.9*lz+0.0001;
                    zpMinus = 0.1*lz-0.0001;
                else
                    zpPlus = 0;
                    zpMinus = lz;
                end                
            end
            
            % check closest distance in z direction               
            if abs(z0-zpPlus) < abs(z0 - zpMinus)
                zp0 = zpPlus;
            else
                zp0 = zpMinus;
            end
            
            % Modify ys: (j-+2,k),(j-+1,k),(j,k-+2),(j,k-+1)         
            deltay=round(abs(y0-yp0),10);
            deltays=y0-yp0;
            deltaz=round(abs(z0-zp0),10);
            deltazs=z0-zp0;
            
            if deltaz<dz
                if round(deltaz,mz+2)~=0
                    zs(2-sign(deltazs))=z0-sign(deltazs)*deltaz;
                end
            end           
            if deltay<dy
                if round(deltay,my+2)~=0 
                    ys(2-sign(deltays))=y0-sign(deltays)*deltay;                     
                end
            end 
            
            bj=zeros(1,3);
            bk=zeros(1,3);          
            bj(1)=2/((ys(2)-ys(1))*(ys(3)-ys(1)));
            bj(2)=-2/((ys(3)-ys(2))*(ys(2)-ys(1))); 
            bj(3)=2/((ys(3)-ys(2))*(ys(3)-ys(1)));
            bk(1)=2/((zs(2)-zs(1))*(zs(3)-zs(1)));
            bk(2)=-2/((zs(3)-zs(2))*(zs(2)-zs(1))); 
            bk(3)=2/((zs(3)-zs(2))*(zs(3)-zs(1)));  
            
            index =(k-1)*ny+j;
            A(index,index)= bj(2) + bk(2);        
            if abs(deltay)<dy           
                if j==1
                  A(index,index+ny-1)= bj(1);
                  A(index,index+1)= 0;
                  
                elseif j==ny
                  A(index,index-ny+1)= bj(3);
                  A(index,index-1)= 0;

                else
                  A(index,index-sign(deltays))=0;       
                  A(index,index+sign(deltays))= bj(2+sign(deltays));

                end
            end                      
            if abs(deltaz)<dz 
                if k==1
                  A(index,index+(nz-1)*ny)= bk(1);
                  A(index,index+ny)= bk(3);
                elseif k==nz
                  A(index,index-ny)= bk(3);
                  A(index,index-(nz-1)*ny)= bk(1);
                else
                  A(index,index-sign(deltazs)*ny)= bk(2-sign(deltazs));
                  A(index,index+sign(deltazs)*ny)= bk(2+sign(deltazs));
                end
            end
        end
    end
end
            
%% No slip at bottom condition
for k = 1:nz
    j = 1;
    index = (k-1)*ny+j ;
    A(index,:) = 0;
    A(index,index) = 1/dy^2;
    P(index)=0;
end

%% Solution
figure
spy(A)
u = mldivide(A,P);
u = reshape(u,ny,nz); % converts u into size (ny,nz)

figure
surfc(z,y,u)
title('u profile')
ylabel('y', 'FontSize', 18)
xlabel('z', 'FontSize', 18)
set(gca,'FontSize',18)

figure
pcolor(u)
shading interp
axis image
colorbar
set(gca,'layer','top')
title('u')
colorbar
ylabel('y')
xlabel('z')

% %% Analytical solution for semi-circlular riblets
% % Used to validate results
% za = linspace(0,lz,nz);
% ya = linspace(0,ly,ny);
% nu = 1;
% R = lz/2;
% index = logical(S);
% for j=1:length(ya)
%     for k=1:length(za)
%         f(j,k)=0.25*(R^2-((ya(j)-ly)^2+(za(k)-lz/2)^2));
% %         fr(j,k)=0.25*(R^2-((ya(j))^2+(za(k))^2));
%     end
% end
% Ubmaen=mean(mean(f));
% % Ubmaenr=mean(mean(fr));
% ua=-1/nu*dpdx*f;
% ua(index)=0;
% umax=max(max(ua));
% umaxmax=max(max(u));
% umsa = mean(ua(ny,:),2);

% %% Verify
% for k = 1:nz
%     Sxtot(k) = (u(ny,k)-u(ny-1,k))/dy;
% end
% Sxtotm = mean(Sxtot);

% %% average in z
% ums = mean(u(ny,:),2);


%% calculate coefficients
integraly = dy*trapz(full(u)); % integrate u*dy with z constant
F1 = dz*trapz(integraly) % integration in y and z

u_surf = u(end,:);
F2 = dz*trapz(full(u_surf)) % integration in z at top layer

if savefile == 1
    tit = ['No_points_' num2str(ny) '_Shape_' num2str(shape) '.mat'];
    save(tit,'ny','relError', 'relErrorint', 'integralz');
end

toc