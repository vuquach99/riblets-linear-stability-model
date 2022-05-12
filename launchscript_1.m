% Linear stability analysis for one riblet spacing over multiple wavelengths
clear
close all

%% Geometry parameters
s = 50; % riblet spacing s+
shape = 'circle';
G1 = 0.008056486729527 % 0.016764297863788 
G2 = 0.072447438454782 % 0.150404461915805
F1 = 0.002829009694317
F2 = 0.016652818813973

%% Other inputs
savefile = 1;
Rt = 550; % friction Reynolds number
% Sweep through a range of wavelengths
nosmod = 256; % number of modes

% wavelength parameters (lxp = friction lambda)
nx = 100; % number of wavelengths
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
most_imag = [];
most_real = [];
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
    fname = ['Rt' num2str(Rt) '_' shape '_Lvp' num2str(Lvpp) '_Ny' num2str(nosmod) '.mat']
    ribstab
    if savefile == 1
        save(fname,'Rt','Max_unstab','y','Lvpp','nosmod','ut','lxp','maxeigvc','maxeigvl','eigvals')
    end
    most_imag(end+1) = imag(Most_unstab)/ut/Rt; % /utau/Rt for channel units
    most_real(end+1) = real(Most_unstab)/ut/Rt;
end

%% Plots most amplified mode of each wavenumber
% figure(1)
% hold on
% for r = 1:size(Lvp,1)
%     Lvpp = Lvp(r);      
%     omega_imag = imag(Max_unstab)/ut/Rt; % /utau/Rt for channel units
%     omega_real = real(Max_unstab)/ut/Rt;
%     plot(omega_real,omega_imag,'.',...
%         'Color', [(r-1)*1/(length(Lvp))', 0, 1-(r-1)*1/(length(Lvp))],'MarkerSize', 10)    
% end
% set(gcf,'position',[160 280 800 600])
% set(gca,'Xlim',[-1 2])
% set(gca,'Ylim',[-1.5 0.5])
% set(gca,'Fontn','Times','FontSize',18,'LineWidth',2)
% xlabel('\omega_r','FontAngle','italic')
% ylabel('\omega_i','FontAngle','italic')
% title('Orr-Sommerfeld spectrum')
% box on