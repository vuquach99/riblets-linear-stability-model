% Linear stability analysis over one wavelength
clear
% close all

%% Geometry parameters
s = [5 10 20 30 40 50]'; % riblet spacing s+
s = 50;
shape = 'testcircle';
G1 = 0.016764297863788 % 0.016764297863788 
G2 = 0.150404461915805 % 0.150404461915805
F1 = 0.012273817101988 % analytical = 0.012271846303085
F2 = 0.041563532761828 % analytical = 0.041666666666667

%% Other inputs
savefile = 0;
Rt = 550; % friction Reynolds number
% Sweep through a range of wavelengths
nosmod = 256; % number of modes

% wavelength parameters (lxp = friction lambda)
nx = 1; 
lxpmin = 9000;
lxpmax = 9001;

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

%% Main loop, plots Orr-Sommerfeld spectrum
omega_imag = zeros(nosmod+1,size(Lvp,1));
omega_real = zeros(nosmod+1,size(Lvp,1));
most_imag = [];
most_real = [];
figure
hold on

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
    omega_imag(:,jK) = imag(eigvals)/ut/Rt; % /ut/Rt for channel units
    omega_real(:,jK) = real(eigvals)/ut/Rt;
    most_imag(end+1) = imag(Most_unstab)/ut/Rt;
    most_real(end+1) = real(Most_unstab)/ut/Rt;
    
    plot(omega_real(:,jK),omega_imag(:,jK),'.',...
        'Color',[(jK-1)*1/(length(Lvp))',0,1-(jK-1)*1/(length(Lvp))],'MarkerSize', 10)
    plot(most_real(jK),most_imag(jK),'d',...
        'Color',[(jK-1)*1/(length(Lvp))',0,1-(jK-1)*1/(length(Lvp))],'MarkerSize', 10)
end
set(gcf,'position',[160 280 800 600])
set(gca,'Xlim',[-0.5 2.5])
set(gca,'Ylim',[-2 0.5])
% set(gca,'Xlim',[0 0.1])
% set(gca,'Ylim',[-0.1 0.05])
set(gca,'Fontn','Times','FontSize',18,'LineWidth',2)
xlabel('\omega_r','FontAngle','italic')
ylabel('\omega_i','FontAngle','italic')
legend('','Most amplified mode','FontSize',18)
title('Orr-Sommerfeld spectrum')
box on