% Plots the change in lg_0.06 when boundary is moved
clear
close all

figure(1)
hold on
% Triangle3:
sigma_max = 0.060562571920025;
lg_06 = [12.7*sqrt(0.5 + sqrt(3)/4) 9.85];
l_U = [0.194168517734770/(1+sqrt(3)/2) 0.194168517734770/(1+sqrt(3)/2)]; 
plot(l_U,lg_06,'--^','Color','#003688','LineWidth',2,'MarkerSize',18)
plot(l_U(1),lg_06(1),'^','Color','#003688','LineWidth',2,'MarkerSize',18,'MarkerFaceColor','#003688')

% Triangle6:
sigma_max = 0.061947584726442;
lg_06 = [21.27*sqrt(sqrt(3))/2 11.6];
l_U = [0.165528802094570/(sqrt(3)/2) 0.165528802094570/(sqrt(3)/2)];
plot(l_U,lg_06,'-->','Color','#9364CD','LineWidth',2,'MarkerSize',18)
plot(l_U(1),lg_06(1),'>','Color','#9364CD','LineWidth',2,'MarkerSize',18,'MarkerFaceColor','#9364CD')

% Triangle9:
sigma_max = 0.063644988097432;
lg_06 = [35.6/2 12.5];
l_U = [0.139389989079762/0.5 0.139389989079762/0.5];
plot(l_U,lg_06,'--v','Color','#E32017','LineWidth',2,'MarkerSize',18)
plot(l_U(1),lg_06(1),'v','Color','#E32017','LineWidth',2,'MarkerSize',18,'MarkerFaceColor','#E32017')

% Semicircle:
sigma_max = 0.065037357706810;
lg_06 = [27.128*sqrt(pi/8) 11.25];
l_U = [0.202240488839837/0.5 0.202240488839837/0.5];
plot(l_U,lg_06,'--o','Color','#E32017','LineWidth',2,'MarkerSize',18)
plot(l_U(1),lg_06(1),'o','Color','#E32017','LineWidth',2,'MarkerSize',18,'MarkerFaceColor','#E32017')

% Trapezium:
sigma_max = 0.059962042273510;
lg_06 = [21.27*sqrt(sqrt(3))/2 10.8];
l_U = [0.165528802094570/0.5 0.165528802094570/0.5];
plot(l_U,lg_06,'--d','Color','#9364CD','LineWidth',2,'MarkerSize',18)
plot(l_U(1),lg_06(1),'d','Color','#9364CD','LineWidth',2,'MarkerSize',18,'MarkerFaceColor','#9364CD')

% Blade:
sigma_max = 0.061195961556218;
lg_06 = [19*sqrt(2/5) 10.8];
l_U = [0.112669343050277/0.5 0.112669343050277/0.5];
plot(l_U,lg_06,'--+','Color','#003688','LineWidth',2,'MarkerSize',18)
plot(l_U(1),lg_06(1),'+','Color','#003688','LineWidth',9,'MarkerSize',18,'MarkerFaceColor','#003688')

set(gcf,'position',[160 280 800 600])
set(gca,'Ylim',[9 19])
set(gca,'Xlim',[0 0.45])
set(gca,'Fontn','Times','FontSize',32,'LineWidth',2)
xlabel('$\ell_U^+/k^+$','Interpreter','latex','FontSize',40)
ylabel('$\ell_{g,0.06}^+$','Interpreter','latex','FontSize',40)
box on