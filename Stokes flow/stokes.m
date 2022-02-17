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
% Sides of the domain
ly = 0.5; 
lz = 1;

% Define grid sizes
ny = 30; % number of points y (wall-normal) domain is enclosed
nz = 60; % number of points in z (spanwise) % periodic in z-direction
dy = ly/(ny-1); % length of sub-intervals in y-axis
dz = lz/(nz-1); % length of sub-intervals in z-axis
z = linspace(0,lz,nz);
y = linspace(0,ly,ny);

dpdx = -1; % Pressure gradient
Sx = 0; % Shear at the top, could set to 1 (if normalised, always 1)

geometry = 2;
% 0 = parabola, 1 = triangle, 2 = semi-circle
% 3 = trapezium, 4 = blades, 5 = ellipse

savefile = 0;

%% Build S matrix - Grid Points of riblets
% S = 1 for points within and on boundary and = 0 elsewhere 
S = zeros(ny,nz);
if geometry == 0 % parabola
    parabola = (2*(z-lz/2)).^2;
    shape = 'parabola';
    for k = 1:nz
        for j = 1:ny
            if y(j)<=parabola(k)
                S(j,k) = 1;
            else
                S(j,k) = 0;
            end
        end
    end
end

if geometry == 1 % triangle
    shape = 'triangle';
    z = linspace(0,1,nz);
    triangle = (z-lz/2);%+1+0.5*z)+0.5*lz;
    triangle2 = -z+lz/2;
    % triangle = (2*z-lz/2);%+1+0.5*z)+0.5*lz;
    % triangle2 = -2*z+lz/2;
    for k = 1:nz
        for j = 1:ny
            if y(j)<=triangle(k) || y(j)<=triangle2(k)
                S(j,k) = 1;
            else
                S(j,k) = 0;
            end
        end
    end
end

if geometry == 2 % semi-circle
    circle = (-sqrt((lz/2)^2-(z-lz/2).^2))+ly;
    shape = 'circle';
    for k = 1:nz
        for j = 1:ny
            if y(j)<=circle(k)
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

if geometry == 3 % trapezium
    shape='trapezium';
    z= linspace(0,1,nz);
    trapezium = (1.333*z-1.6667*lz/2);%+1+0.5*z)+0.5*lz;
    trapezium2 = (-1.333*z+lz/2);
    for k=1:nz
        for j=1:ny
            if y(j)<=trapezium(k) || y(j)<=trapezium2(k)
                S(j, k) = 1;
            else
                S(j, k) = 0;
            end
        end
    end
    S(1,:) = ones(1,nz);
end

%% Build Sd matrix
% Sd = 1 for points within and on boundary and = 0 elsewhere 
% riblets now bigger by d = max(dy,dz)
Sd = zeros(ny,nz);
if geometry == 0 % parabola
    z = linspace(0,1,nz);
    parabola_d = (2*(z-lz/2)).^2+dz;
    for k = 1:nz
        for j = 1:ny
            if y(j)<parabola_d(k)
                Sd(j,k) = 1;
            else
                Sd(j,k) = 0;
            end
        end
    end
end

if geometry == 1 % triangle
    z= linspace(0,1,nz);
    triangle = (z-lz/2)+dy;%+1+0.5*z)+0.5*lz;
    triangle2 = -z+lz/2+dy;
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
    z = linspace(0,1,nz);
    circle = (-sqrt((lz/2)^2-(z-lz/2).^2))+ly+dy;
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
    z= linspace(0,1,nz);
    trapezium = (1.333*z-1.6667*lz/2+dy);%+1+0.5*z)+0.5*lz;
    trapezium2 = (-1.333*z+lz/2+dy);
    for k=1:nz
        for j=1:ny
            if y(j)<trapezium(k) || y(j)<trapezium2(k)
                Sd(j, k) = 1;
            else
                Sd(j, k) = 0;
            end
        end
    end
    Sd(1:2,:) = ones(2,nz);
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
            z0 = k*dz ;          
            my = abs(floor(log10(dy)));
            mz = abs(floor(log10(dz)));
            
        %ys: (j-+2,k);(j-+1,k);(j,k-+2);(j,k-+1) 
            ys=[y0-dy, y0, y0+dy];
            zs=[z0-dz, z0, z0+dz];  
            
            if geometry == 0 % parabola
                yp0 = (2*(z0-lz/2)).^2;
                zpPlus = sqrt(y0)/2+lz/2;
                zpMinus = -zpPlus;
            end
            
            if geometry == 1 % triangle
                if z0 <= lz/2
                    yp0 = z0 - lz/2;
                else
                    yp0 = -z0+lz/2;
                end
                zpPlus = lz/2 - y0;
                zpMinus = y0 + lz/2;
            end
            
            if geometry == 2 % semi-circle
                yp0=-(-sqrt((lz/2)^2-(z0-lz/2).^2))+ly;
                zpPlus=sqrt((lz/2)^2-(y0-ly).^2)+lz/2;
                zpMinus=-zpPlus;
            end
            
            if geometry == 3 % trapezium
                if z0 <= lz/2
                    yp0 =  (-1.333*z+lz/2);
                else
                    yp0 = (1.333*z-1.6667*lz/2);
                end
                zpPlus = (lz/2 - y0)/1.333;
                zpMinus = (1.6667*lz/2 - y0)/1.333;
                yp0 = -yp0;
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
u = reshape(mldivide(A,P),ny,nz);
Ub = mean(mean(u));
u_max = full(max(u(:)));
if dpdx == -1 && Sx == 0 && geometry == 2
    fprintf('4*u_max = %f\n', u_max*4)
else
    fprintf('u_max = %f\n', u_max)
end
% Lws=Ub/s;
figure
surfc(z,y,u)
% daspect([1 1 1])
% axis tight
% axis fill
% axis image
ylabel('y', 'FontSize', 18)
xlabel('z', 'FontSize', 18)
set(gca,'FontSize',18)
% us = [zeros(1,nz); u];
us = [u(:,end) u];

%% Analytical solution for semi-circlular riblets
% Used to validate results
za = linspace(0,lz,nz);
ya = linspace(0,ly,ny);
nu = 1;
R = lz/2;
index = logical(S);
for j=1:length(ya)
    for k=1:length(za)
        f(j,k)=0.25*(R^2-((ya(j)-ly)^2+(za(k)-lz/2)^2));
%         fr(j,k)=0.25*(R^2-((ya(j))^2+(za(k))^2));
    end
end
Ubmaen=mean(mean(f));
% Ubmaenr=mean(mean(fr));
ua=-1/nu*dpdx*f;
ua(index)=0;
umax=max(max(ua));
umaxmax=max(max(u));

%% Verify
for k = 1:nz
    Sxtot(k) = (u(ny,k)-u(ny-1,k))/dy;
end
Sxtotm = mean(Sxtot);

%% average in z
ums = mean(u(ny,:),2);
umsa = mean(ua(ny,:),2);

%% Error and Grid Independence
fint=f;
fint(index)=0;
integraly = dy*trapz(full(u));
integralz = dz*trapz(integraly) % integration in y and z
% f_integraly=dy*trapz(full(fint));
% f_integralz=dz*trapz(f_integraly)
% relError=abs((ums-umsa)/umax*100);
% relErrorint=abs((0.0122602323037932-integralz)/umax*100);
relErrorint = abs((0.0122602323037932-integralz)/0.0122602323037932*100);
relError = abs((ua-u)/umax*100);
points_y = ny;

if savefile == 1
    tit = ['No_points_' num2str(ny) '_Shape_' num2str(shape) '.mat'];
    save(tit,'points_y','relError', 'relErrorint', 'integralz');
end

toc