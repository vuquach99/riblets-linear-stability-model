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
a = 1.9; % aspect ratio ly:lz
lz = 1;
ly = a*lz;

% Define grid sizes
n = 25;
nz = 2*n; % number of points in z (spanwise) % periodic in z-direction
ny = a*nz; % number of points y (wall-normal) domain is enclosed
dy = ly/(ny-1); % length of sub-intervals in y-axis
dz = lz/nz; % length of sub-intervals in z-axis
z = (0:nz-1)*dz;
y = linspace(0,ly,ny);

dpdx = 0; % Pressure gradient
Sx = 1; % Shear at the top, could set to 1 (if normalised, always 1)

geometry = 1; 
% 0 = parabola 
% 1 = triangle (k/s = 1.866 for 30deg, 0.5 for 90deg, 0.866 for 60deg)
% 2 = semi-circle (k/s = 0.5)

angle = 90; % 30/60/90 degrees, for triangles only
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
    if angle == 90        
        shape = 'triangle9';
        triangle = z-lz/2;
        triangle2 = -z+lz/2;
    end
    if angle == 60
        shape = 'triangle6';
        triangle = sqrt(3)*z-sqrt(3)*lz/2;
        triangle2 = -sqrt(3)*z+sqrt(3)*lz/2;
    end
    if angle == 30
        shape = 'triangle3';
        triangle = (2+sqrt(3))*z-1-sqrt(3)*lz/2;
        triangle2 = -(2+sqrt(3))*z+1+sqrt(3)*lz/2;
    end
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
figure
spy(S)
title('S matrix')

%% Build Sd matrix
% Sd = 1 for points within and on boundary and = 0 elsewhere 
% riblets now bigger by d = max(dy,dz)
Sd = zeros(ny,nz);
if geometry == 0 % parabola
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
    if angle == 90
        triangle = z-lz/2+dy;
        triangle2 = -z+lz/2+dy;
    end
    if angle == 60
        triangle = sqrt(3)*z-sqrt(3)*lz/2+dy;
        triangle2 = -sqrt(3)*z+sqrt(3)*lz/2+dy;
    end
    if angle == 30
        triangle = (2+sqrt(3))*z-1-sqrt(3)*lz/2+dy;
        triangle2 = -(2+sqrt(3))*z+1+sqrt(3)*lz/2+dy;
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
    circle = (-sqrt((lz/2)^2-(z-lz/2).^2))+lz/2+dy;
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
figure
spy(Sd)
title('Sd matrix')

Sc = Sd-S;
figure
spy(Sc)
title('Sc matrix')