% Main linear stability analysis programme.
clear
close all

%% Geometry parameters
s = 50; % riblet spacing
% s = [1 2 5 10 11 14 16 18 20 25 30 40 50]'; %[5 10 20 30 40 50]';

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
Lwp = F1^(1/3)*s;
Lsp = F2^(1/2)*s;
Lhq = G1^(1/2)*s;
Lsq = G2*s;

[D0,D1,D2,D3,D4] = dmat(nosmod); % Chebyshev polynomials

% Cess turbulent velocity profile inside channel
[y,nut,uc,Re] = turprof_generic(nosmod,Aint,kapa,eddyfrac,Rt); 
ut = Rt/Re; % friction velocity

% wavelength vector
lxp = log(lxpmin):(log(lxpmax)-log(lxpmin))/(nx-1):log(lxpmax);
lxp = fliplr(exp(lxp));
alp0 = 2*pi*Rt./lxp; % wavenumber in channel units

%% Main loop
for jK = 1:size(Lwp,1)
    Lwpp=Lwp(jK); % wall-normal coeff. in wall units
    Kw=(Lwpp/Rt)^3; % wall normal coeff. in outer units cubed (Lw3)
    iKw=1/Kw/Re;

    Lspp=Lsp(jK); % streamwise coeff. in wall units
    Sw=(Lspp/Rt)^2; % streamwise coeff. in outer units squared (Ls2)
    iSw=1/Sw/Re;

    % Coefficients for shear driven flow
    lhq=(Lhq(jK)/Rt)^2; 
    lsq=Lsq(jK)/Rt;
    L=lhq/Kw/Re;

    ribstab
    fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat'];
    save(fname,'Rt','Max_unstab','y','Lwpp','nosmod','ut','lxp','maxeigvc','maxeigvl','eigvals')
end

%% Plotting
n = length(Lwp);
% Orr-Sommerfeld spectrum
figure(1)
for r=1:length(Lwp)
    Lwpp=Lwp(r);
    fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat']
    load(fname)
    a{r}= imag(eigvals)/ut/Rt>0.2;
    b{r} = imag(eigvals(a{r}));
%     maxim=max(imag(eigvals));
%     maxim=max(maxim)/ut/Rt;
%     plot(real(eigvals)/ut/Rt,imag(eigvals)/ut/Rt,'.', 'Color', [(r-1)*1/n', 0, 1-(r-1)*1/n],'MarkerSize', 10)
     plot(real(eigvals)/ut/Rt,imag(eigvals)/ut/Rt,'.', 'Color', [(r-1)*1/n', 0, 1-(r-1)*1/n],'MarkerSize', 10)
    hold on
end

xlabel('\omega^+_R')
ylabel('\omega^+_I')
title('Orr-Sommerfeld spectrum')
box on
set(gca,'Fontn','Times','FontSize',14,'LineWidth',2)
