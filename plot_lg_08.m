% Plots lg_0.08 when boundary is at l_U
clear
close all

figure(1)
hold on
% Triangle3:
lg_08 = 14.1;
l_U = 0.194168517734770/(1+sqrt(3)/2);
plot(l_U,lg_08,'^','Color','#003688','LineWidth',2,'MarkerSize',18)

% Triangle6:
lg_08 = 16.84;
l_U = 0.165528802094570/(sqrt(3)/2);
plot(l_U,lg_08,'>','Color','#9364CD','LineWidth',2,'MarkerSize',18)

% Triangle9:
lg_08 = 20.4;
l_U = 0.139389989079762/0.5;
plot(l_U,lg_08,'v','Color','#E32017','LineWidth',2,'MarkerSize',18)

% Semicircle:
lg_08 = 19.55;
l_U = 0.202240488839837/0.5;
plot(l_U,lg_08,'o','Color','#E32017','LineWidth',2,'MarkerSize',18)

% Trapezium:
lg_08 = 17.19;
l_U = 0.165528802094570/0.5;
plot(l_U,lg_08,'d','Color','#9364CD','LineWidth',2,'MarkerSize',18)

% Blade:
lg_08 = 14.5;
l_U = 0.112669343050277/0.5;
plot(l_U,lg_08,'+','Color','#003688','LineWidth',2,'MarkerSize',18)

set(gcf,'position',[160 280 800 600])
set(gca,'Ylim',[10 22])
set(gca,'Xlim',[0 0.45])
set(gca,'Fontn','Times','FontSize',32,'LineWidth',2)
xlabel('$\ell_U^+/k^+$','Interpreter','latex','FontSize',40)
ylabel('$\ell_{g,0.08}^+$','Interpreter','latex','FontSize',40)
box on