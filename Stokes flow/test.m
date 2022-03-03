%% Programme to solve for Stokes velocity profile in riblets
% Solve A*u = P
% S = riblets grid points
% P = pressure gradient
% A = Laplacian coefficients
% Sc = Sd - S = riblets surface grid points


close all
tic

%% Input data
% Sides of the periodic-in-z domain
ly = 1.5; % aspect ratio ly:lz
lz = 1; % always

% Define grid sizes
n = 40;
nz = 2*n; % number of points in periodic z (spanwise)
ny = ly*nz; % number of points y (wall-normal) domain is enclosed
dy = ly/(ny-1); % length of sub-intervals in y-axis
dz = lz/nz; % length of sub-intervals in z-axis
z = (0:nz-1)*dz;
y = linspace(0,ly,ny);

dpdx = 0; % Pressure gradient
Sx = 1; % Shear at the top, could set to 1 (if normalised, always 1)

geometry = 4; 
% 0 = parabola (k/s = 1)
% 1 = triangle (k/s = 1.866 for 30deg, 0.5 for 90deg, 0.866 for 60deg)
% 2 = semi-circle (k/s = 0.5)
% 3 = trapezium (k/s = 0.5; tip half-angle = 15deg)
% 4 = blade (k/s = 0.5; t/s = 0.2)
angle = 90; % 30/60/90 degrees, for triangles only

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

if geometry == 4 % blade
    shape = 'blade';
    height = 0.5;
    for j = 1:height*nz
        for k = 1:nz
            if z(k) <= 0.1*lz || z(k) >= 0.9*lz
                S(j,k) = 1;
            else
                S(j,k) = 0;
            end
        end
    end
end

% figure
% spy(S)
% title('S matrix')

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

if geometry == 4 % blade
    for j = 1:height*nz+1
        for k = 1:nz
            if z(k) <= 0.1*lz+dz || z(k) >= 0.9*lz-dz
                Sd(j,k) = 1;
            else
                Sd(j,k) = 0;
            end
        end
    end
    Sd(1,:) = ones(1,nz);
end

% figure
% spy(Sd)
% title('Sd matrix')

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
            y0 = (j-ny)*dy; %?
            z0 = k*dz ;          
            my = abs(floor(log10(dy)));
            mz = abs(floor(log10(dz)));            
            % ys: (j-+1,k),(j,k-+1) 
            ys=[y0-dy, y0, y0+dy];
            zs=[z0-dz, z0, z0+dz];
            
            % find the points of intersection
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
                    zpPlus = y0 + lz/2;
                    zpMinus = -y0 + lz/2;
                end
            end
            
            if geometry == 4 % blade
                if z0 <= 0.1*lz || z0 >= 0.9*lz
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
            if abs(z0-zpPlus) < abs(z0-zpMinus)
                zp0 = zpPlus;
            else
                zp0 = zpMinus;
            end
            
            % Modify ys: (j-+1,k),(j,k-+1)         
            deltay=round(abs(y0-yp0),10);
            deltays=y0-yp0;
            deltaz=round(abs(z0-zp0),10);
            deltazs=z0-zp0;            
            if deltaz<dz % z'
                if round(deltaz,mz+2)~=0
                    zs(2-sign(deltazs))=z0-sign(deltazs)*deltaz;
                end
            end           
            if deltay<dy % y'
                if round(deltay,my+2)~=0 
                    ys(2-sign(deltays))=y0-sign(deltays)*deltay;                     
                end
            end 
            
            % finite difference
            bj = zeros(1,3);
            bk = zeros(1,3);
            
            bj(1) = 2/((ys(2)-ys(1))*(ys(3)-ys(1)));
            bj(2) = -2/((ys(3)-ys(2))*(ys(2)-ys(1))); 
            bj(3) = 2/((ys(3)-ys(2))*(ys(3)-ys(1)));
            
            bk(1) = 2/((zs(2)-zs(1))*(zs(3)-zs(1)));
            bk(2) = -2/((zs(3)-zs(2))*(zs(2)-zs(1))); 
            bk(3) = 2/((zs(3)-zs(2))*(zs(3)-zs(1)));  
            
            index = (k-1)*ny+j;
            A(index,index) = bj(2) + bk(2);            
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
                    A(index,index+(nz-1)*ny) = bk(1);
                    A(index,index+ny) = bk(3);
                elseif k==nz
                    A(index,index-ny)= bk(3);
                    A(index,index-(nz-1)*ny)= bk(1);
                else
                    A(index,index-sign(deltazs)*ny)= bk(2-sign(deltazs)); %?
                    A(index,index+sign(deltazs)*ny)= bk(2+sign(deltazs));
                end
            end
        end
    end
end
A = full(A);            
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
u = reshape(u,ny,nz);
if size(u,1) > height*nz
    u_g = u(1:round(height*nz),:); % velocity profile in groove
    n1 = round(height*nz); % index for plotting
else
    u_g = u;
    n1 = size(u,1);
end
u_max = full(max(u_g(:)));
fprintf('u_max = %f\n', u_max)
figure
surfc(z,y,u)
title('u profile')
axis image
ylabel('y', 'FontSize', 18)
xlabel('z', 'FontSize', 18)
set(gca,'FontSize',18)

figure
surfc(z,y(1:n1),u_g)
title('u profile in groove')

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

%% Verify
for k = 1:nz
    Sxtot(k) = (u(ny,k)-u(ny-1,k))/dy;
end
Sxtotm = mean(Sxtot);

%% average in z and find point of min gradient
ums = mean(u,2);
ums = full(ums);
% find point of min gradient
grad = zeros(1,ny);
for i = 1:ny-1
    grad(i) = (y(i+1)-y(i))/(ums(i+1)-ums(i));
end
grad(ny) = grad(ny-1);
[grad_min, minimum] = min(grad);
boundary = y(minimum)-grad_min*ums(minimum);
y_line = grad_min*ums + boundary;
fprintf('Virtual boundary at y = %f\n', boundary)
% 
% figure
% hold on
% plot(ums, y,'k','Linewidth',1.25)
% plot(ums,y_line,'r','Linewidth',1.25)
% ylabel('$y$','Interpreter','Latex')
% xlabel('$\bar{u}$','Interpreter','Latex')
% set(gca,'FontName','times','FontSize',15)

%% Error and Grid Independence
integraly = dy*trapz(full(u));
integralz = dz*trapz(integraly) % integration in y and z
points_y = ny-1;

if savefile == 1
    tit = ['No_points_' num2str(ny-1) '_Shape_' num2str(shape) '.mat'];
    save(tit,'points_y','relError', 'relErrorint', 'integralz');
end

toc