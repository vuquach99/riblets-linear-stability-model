% Linear stability analysis over multiple wavelengths
tic
clear
close all

%% Geometry parameters
s = 21.27; % riblet spacing s+
% s = [5 7.5 10 12.5 15 17.5 20 30 40 50]';

shape = 'trapezium1_s';
% at virtual origin
G1 = 0.009516234685273
G2 = 0.072744932298483
F1 = 0.004823931048779
F2 = 0.022935819820992

% F1 = 0.002829009694317
% F2 = 0.016652818813973

%% Other inputs
savefile = 0;
Rt = 550; % friction Reynolds number
% Sweep through a range of wavelengths
nosmod = 512; % number of modes

% wavelength parameters (lxp = friction lambda)
nx = 200; % number of wavelengths
lxpmin = 10;
lxpmax = 10000;

% Parameters for generating velocity profile
kapa = 0.426;
Aint = 25.4;
eddyfrac = 1;

%% Calculates essential parameters
% permeability expressions
Lvp = F1^(1/3)*s; % Lwp
Lup = F2^(1/2)*s; % Lsp
Lvs = G1^(1/2)*s; % Lhq
Lus = G2*s; % Lsq

[D0,D1,D2,D3,D4] = dmat(nosmod); % Chebyshev polynomials

% Cess turbulent velocity profile inside channel
[y,nut,U,Re] = turprof_generic(nosmod,Aint,kapa,eddyfrac,Rt); 
ut = Rt/Re; % friction velocity

% wavelength vector
lxp = log(lxpmin):(log(lxpmax)-log(lxpmin))/(nx-1):log(lxpmax);
lxp = fliplr(exp(lxp));
alp0 = 2*pi*Rt./lxp; % wavenumber in channel units

%% Main loop
for jK = 1:size(s,1)
    sp = s(jK);
    % Pressure-driven coefficients    
    Lvpp = Lvp(jK);
    Kvp = (Lvpp/Rt)^3; % wall normal coefficients in outer units cubed
    Lupp = Lup(jK);
    Kup = (Lupp/Rt)^2; % streamwise coefficients in outer units squared
    % Shear-driven coefficients
    Kvs = (Lvs(jK)/Rt)^2; % Lvs
    Kus = Lus/Rt; % Lus
    % Calculates stuff
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat']
    ribstab
    if savefile == 1
        save(fname,'Rt','Most_unstab','Most_lxp','Max_unstab','y','sp','nosmod',...
            'ut','lxp','maxeigvc','maxeigvl','eigvals')
    end
end
toc