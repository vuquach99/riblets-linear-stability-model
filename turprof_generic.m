function [y,nut,u,ReB] = turprof_generic(N,Aint,kapa,eddyfrac,Rt);
%
% Function to create velocity and eddy viscosity profiles to pass to turchan using Chebyshev
% pseudospectral discretization for plane turbulent channels.
%
% INPUT
%
%
% OUTPUT
% 
% A,B		= Orr-Sommerfeld matrices   ==>   d = inv(B)*A;
%

% mean velocity and turbulent viscosity
% uses analytical expression from Reynolds and Tiedermann
%
%%%%%%%%%%%%%  use full eddy viscosity to compute the profile 
y = cos(pi*(0:N)/N)';  Nos=N+1;
%the cess eddy viscosity model
y0 = 1-y;
nut=0.5*(1 + kapa^2*Rt^2/9*(2*y0-y0.^2).^2.*(3-4*y0+2*y0.^2).^2.*...
    (1-exp(-y0*Rt/Aint)).^2).^.5+0.5;
%calculamos el perfil medio
u=Rt*cumtrapz(y0,(1-y0)./nut); u=u-u(1);
du  = -y./nut;
%%%     simetriza
n0  = fix(Nos/2);
n1  = Nos-n0;
du  = [du(1:n1); -du(n0:-1:1)];
nut = [nut(1:n1); nut(n0:-1:1)];
u   = [u(1:n1); u(n0:-1:1)];
%%%     reescala
Ub  = trapz(y0,u)/trapz(y0,ones(Nos,1));
u   = u/Ub;
ReB = Ub*Rt
nut = nut/ReB;
du  = du*Rt/Ub;

upp=du*0;
upp(2:end-1) = (du(3:end)-du(1:end-2))./(y(3:end)-y(1:end-2));
upp(1)   = (du(2)-du(1))/(y(2)-y(1));
upp(end) = upp(1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  decrease now the eddy viscosity by a factor eddyfrac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nut = nut(1)+(nut-nut(1))*eddyfrac;%for cess formula

% nut = 1+nut*eddyfrac;%for the "JC" numerical model

%the not viscosity in the center model
%nut=cosh(evalpha*y)./cosh(evalpha);%alphaev is the shape parameter.
%
%plot(y,nut)
%%%%%%%%%%%% read of the nut "JC" numerical expresion
% fid=fopen('../nutjctch.bin','r');
% a=fread(fid,1,'int');
% nut=fread(fid,Nos,'float64');%y en el mallado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nutp=nut*0;
nutp(2:end-1) = (nut(3:end)-nut(1:end-2))./(y(3:end)-y(1:end-2));
nutp(1)   = (nut(2)-nut(1))/(y(2)-y(1)); 
nutp(end) = -nutp(1);
%
nutpp=nut*0;
nutpp(2:end-1) = (nutp(3:end)-nutp(1:end-2))./(y(3:end)-y(1:end-2));
nutpp(1)   = (nutp(2)-nutp(1))/(y(2)-y(1));
nutpp(end) = nutpp(1);
%
nutppp=nut*0;
nutppp(2:end-1) = (nutpp(3:end)-nutpp(1:end-2))./(y(3:end)-y(1:end-2));
nutppp(1)   = (nutpp(2)-nutpp(1))/(y(2)-y(1));
nutppp(end) = -nutppp(1);

%%%%%%%%%%%%%  collect into something manageable %%%%%%%%
u=[u,du,upp];
nut=[nut,nutp,nutpp,nutppp];
