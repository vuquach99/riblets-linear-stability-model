% Main linear stability analysis programme.
clear
close all

%% Geometry parameters
s = 50; % riblet spacing
% s = [1 2 5 10 11 14 16 18 20 25 30 40 50]'; %[5 10 20 30 40 50]';
shape = 'circle';
G1 = 0; % 0.041857060995110
G2 = 0; % 0.234528106122330
F1 = 0.012273817101988; % analytical = 0.012271846303085
F2 = 0.041563532761828; % analytical = 0.041666666666667

%% Other inputs
Rt = 550; % friction Reynolds number
% Sweep through a range of wavelengths
nosmod = 256; % number of modes

% wavelength parameters (lxp = friction lambda)
nx = 100; % number of wavelengths
lxpmin = 10;
lxpmax = 1000;
% for investigating a single wavelength value
% nx = 1; 
% lxpmin=56;
% lxpmax=57;

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
for jK = 1:size(Lvp,1)
    % Pressure-driven coefficients
    Lvpp = Lvp(jK);
    Kvp = (Lvpp/Rt)^3; % wall normal coefficients in outer units cubed
    Lupp = Lup(jK);
    Kup = (Lupp/Rt)^2; % streamwise coefficients in outer units squared
    % Shear-driven coefficients
    Kvs = (Lvs(jK)/Rt)^2; % Lvs
    Kus = Lus/Rt; % Lus
    % Calculates stuff
    fname = ['ribstab_Rt' num2str(Rt) '_' shape '_Ny' num2str(nosmod) '.mat'];
    ribstab
    save(fname,'Rt','Max_unstab','y','Lvpp','nosmod','ut','lxp','maxeigvc','maxeigvl','eigvals')
end

%% Plots Orr-Sommerfeld spectrum
figure(1)
for r = 1:length(Lvp)
    Lvpp = Lvp(r);
    fname = ['ribstab_Rt' num2str(Rt) '_' shape '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    a{r}= imag(eigvals)/ut/Rt>0.2;
    b{r} = imag(eigvals(a{r}));
    plot(real(eigvals)/ut/Rt,imag(eigvals)/ut/Rt,'.',...
        'Color', [(r-1)*1/(length(Lvp))', 0, 1-(r-1)*1/(length(Lvp))],'MarkerSize', 10)
    hold on
end

xlabel('\omega^+_R')
ylabel('\omega^+_I')
title('Orr-Sommerfeld spectrum')
box on
set(gca,'Fontn','Times','FontSize',14,'LineWidth',2)