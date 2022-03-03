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
% Sides of the periodic-in-z domain
ly = 0.6; % aspect ratio ly:lz
lz = 1; % always

% Define grid sizes
n = 25;
nz = 2*n; % number of points in periodic z (spanwise)
ny = ly*nz; % number of points y (wall-normal) domain is enclosed
dy = ly/(ny-1); % length of sub-intervals in y-axis
dz = lz/nz; % length of sub-intervals in z-axis
z = (0:nz-1)*dz;
y = linspace(0,ly,ny);



dpdx = -1; % Pressure gradient
Sx = 0; % Shear at the top, could set to 1 (if normalised, always 1)

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
Sc = Sd-S;
figure
spy(Sc)
title('Sc matrix')