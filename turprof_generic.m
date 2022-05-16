function [y,nut,U,ReB] = turprof_generic(N,Aint,kapa,eddyfrac,Rt)
% Creates velocity and eddy viscosity profiles to pass to turchan.
% Uses Chebyshev pseudospectral discretization for plane turbulent channels.
% ' = d/dy

% INPUT:
% N = number of modes
% Aint, kapa, eddyfrac = parameters
% Rt = friction Reynolds number

% OUTPUT:
% y = wall-normal coordinate
% U = mean velocity profile
% nut = full eddy viscosity to compute the profile
% Rt = bulk Reynolds number

% uses analytical expression from Reynolds and Tiedermann
y = cos(pi*(0:N)/N)';
Nos = N+1;

% Cess eddy viscosity model
y0 = 1-y;
% y0 = y0 + (spanwise protrusion height, Luchini in outer units)
% remember to set nut to 1 below y0
nut = 0.5*(1 + kapa^2*Rt^2/9*(2*y0-y0.^2).^2.*(3-4*y0+2*y0.^2).^2.*...
    (1-exp(-y0*Rt/Aint)).^2).^.5+0.5;

% calculates mean profile
u = Rt*cumtrapz(y0,(1-y0)./nut);
u = u-u(1);
du = -y./nut;

% symmetrises
n0 = fix(Nos/2);
n1 = Nos-n0;
du = [du(1:n1); -du(n0:-1:1)];
nut = [nut(1:n1); nut(n0:-1:1)];
u = [u(1:n1); u(n0:-1:1)];

% rescales
Ub = trapz(y0,u)/trapz(y0,ones(Nos,1));
u = u/Ub; % U
ReB = Ub*Rt
nut = nut/ReB;
du = du*Rt/Ub; % U'

upp = du*0; % U''
upp(2:end-1) = (du(3:end)-du(1:end-2))./(y(3:end)-y(1:end-2));
upp(1) = (du(2)-du(1))/(y(2)-y(1));
upp(end) = upp(1);

% decreases the eddy viscosity by a factor eddyfrac
nut = nut(1)+(nut-nut(1))*eddyfrac; % nu for Cess' formula

nutp = nut*0; % nu'
nutp(2:end-1) = (nut(3:end)-nut(1:end-2))./(y(3:end)-y(1:end-2));
nutp(1) = (nut(2)-nut(1))/(y(2)-y(1)); 
nutp(end) = -nutp(1);

nutpp = nut*0; % nu''
nutpp(2:end-1) = (nutp(3:end)-nutp(1:end-2))./(y(3:end)-y(1:end-2));
nutpp(1) = (nutp(2)-nutp(1))/(y(2)-y(1));
nutpp(end) = nutpp(1);

nutppp = nut*0; % nu'''
nutppp(2:end-1) = (nutpp(3:end)-nutpp(1:end-2))./(y(3:end)-y(1:end-2));
nutppp(1) = (nutpp(2)-nutpp(1))/(y(2)-y(1));
nutppp(end) = -nutppp(1);

% collects into something manageable
U = [u,du,upp];
nut = [nut,nutp,nutpp,nutppp];
