% Plots lg_0.06 when boundary is at l_U
clear
close all

figure(1)
hold on
% Triangle3:
sigma_max = 0.060562571920025;
lg_06 = 12.7*sqrt(0.5 + sqrt(3)/4);
l_U = 0.194168517734770/(1+sqrt(3)/2);
plot(l_U,lg_06,'^','Color','#003688','LineWidth',2,'MarkerSize',18)

% Triangle6:
sigma_max = 0.061947584726442;
lg_06 = 21.27*sqrt(sqrt(3))/2;
l_U = 0.165528802094570/(sqrt(3)/2);
plot(l_U,lg_06,'>','Color','#9364CD','LineWidth',2,'MarkerSize',18)

% Triangle9:
sigma_max = 0.063644988097432;
lg_06 = 35.6/2;
l_U = 0.139389989079762/0.5;
plot(l_U,lg_06,'v','Color','#E32017','LineWidth',2,'MarkerSize',18)

% Semicircle:
sigma_max = 0.065037357706810;
lg_06 = 27.128*sqrt(pi/8);
l_U = 0.202240488839837/0.5;
plot(l_U,lg_06,'o','Color','#E32017','LineWidth',2,'MarkerSize',18)

% Trapezium:
sigma_max = 0.059962042273510;
lg_06 = 21.27*sqrt(sqrt(3))/2;
l_U = 0.165528802094570/0.5;
plot(l_U,lg_06,'d','Color','#9364CD','LineWidth',2,'MarkerSize',18)

% Blade:
sigma_max = 0.061195961556218;
lg_06 = 19*sqrt(2/5);
l_U = 0.112669343050277/0.5;
plot(l_U,lg_06,'+','Color','#003688','LineWidth',2,'MarkerSize',18)

set(gcf,'position',[160 280 800 600])
set(gca,'Ylim',[10 20])
set(gca,'Xlim',[0 0.45])
set(gca,'Fontn','Times','FontSize',32,'LineWidth',2)
xlabel('$\ell_U^+/k^+$','Interpreter','latex','FontSize',40)
ylabel('$\ell_{g,0.06}^+$','Interpreter','latex','FontSize',40)
box on