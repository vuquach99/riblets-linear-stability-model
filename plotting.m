% Plotting script
clear
close all

%% Maximum amplification growth rate vs wavelength+
figure(1)
hold on

% Set at tips
fname = 'Rt550_circle_Lvp11.5336_Ny256.mat';
load(fname)
nx = length(lxp);
for p = 1:nx
    imag_eigval(p) = imag(Max_unstab(p))/ut/Rt;
end
plot(lxp,imag_eigval,'LineWidth',2)

% Set inside grooves
fname = 'Rt550_circle_Lvp7.0716_Ny256.mat';
load(fname)
nx = length(lxp);
for p = 1:nx
    imag_eigval(p) = imag(Max_unstab(p))/ut/Rt;
end

plot(lxp,imag_eigval,'LineWidth',2)
set(gca,'xscale','log')
yline(0,'--','LineWidth',2)
set(gcf,'position',[160 280 800 600])
set(gca,'Xlim',[10 10000])
set(gca,'Ylim',[-0.3 0.35])
set(gca,'Fontn','Times','FontSize',18,'LineWidth',2)
xlabel('\lambda_x^+','FontAngle','italic');
ylabel('\sigma_{Im}^+','FontAngle','italic');
title('Amplification vs Wavelength for s^+ = 20, semi-circle, no shear')
legend('Set at tips','Set inside grooves','location','Southeast','FontSize',18)
box on