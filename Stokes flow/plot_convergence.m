% Plots convergence results
clear
close all

load('convergence_pressure.mat')

figure(1)
set(gcf,'position',[160 280 800 600])
loglog(n_hist,abs((F1_hist-0.012271846303085)/0.012271846303085)...
    ,'-+','MarkerSize',14,'MarkerIndices',5:10:length(n_hist),'Color','#003688','LineWidth',4)
set(gca,'Fontn','Times','FontSize',32,'LineWidth',4)
ylabel('$err(L^3_{vp})$','Interpreter','latex','FontSize',44)
xlabel('$n_y,n_z$','Interpreter','latex','FontSize',44)
box on

figure(2)
set(gcf,'position',[160 280 800 600])
loglog(n_hist,abs((F2_hist-0.041666666666667)/0.041666666666667)...
    ,'-+','MarkerSize',14,'MarkerIndices',5:10:length(n_hist),'Color','#003688','LineWidth',4)
set(gca,'Fontn','Times','FontSize',32,'LineWidth',4)
ylabel('$err(L^2_{up})$','Interpreter','latex','FontSize',44)
xlabel('$n_y,n_z$','Interpreter','latex','FontSize',44)
box on