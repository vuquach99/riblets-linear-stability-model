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
ly = 0.5; % aspect ratio ly:lz
lz = 1; % always

% Define grid sizes
n = 5;
nz = 2*n+1; % number of points in non-periodic z (spanwise)
ny = ly*(nz-1); % number of points y (wall-normal) domain is enclosed
dy = ly/(ny-1); % length of sub-intervals in y-axis
dz = lz/(nz-1); % length of sub-intervals in z-axis
z = linspace(0,lz,nz);
y = linspace(0,ly,ny);

dpdx = -1; % Pressure gradient
Sx = 0; % Shear at the top, could set to 1 (if normalised, always 1)

geometry = 2;
% 0 = parabola, 1 = triangle, 2 = semi-circle
% 3 = trapezium, 4 = blades

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
    circle = (-sqrt((lz/2)^2-(z-lz/2).^2))+lz/2;
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

if geometry == 4 % blade
    shape = 'blade';
    height = 0.5;
    for j = 1:height*(nz-1)
        for k = 1:nz
            if z(k) < 0.2*lz
                S(j, k) = 1;
            else
                S(j, k) = 0;
            end
        end
        S(j,nz) = 1;
    end
end

figure
spy(S)
title('S matrix')

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
if size(u,1) > n
    u_g = u(1:n,:); % velocity profile in groove
    n1 = n; % index for plotting
else
    u_g = u;
    n1 = size(u,1);
end
u_max = full(max(u(:)));
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