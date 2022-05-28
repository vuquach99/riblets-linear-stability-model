% Plots values of u_g for all riblet shapes
clear
close all

figure(1)
hold on

load('shape_triangle3.mat')
y1 = (y-height)/height;
k = y1 >= 0;
k = find(k,1,'first');
u1 = ums(k);
plot(ums/u1,(y-height)/height,'-^','Color','#003688','Linewidth',2,'MarkerIndices',1:3:length(y1),'MarkerSize',9)

load('shape_triangle6.mat')
y1 = (y-height)/height;
k = y1 >= 0;
k = find(k,1,'first');
u1 = ums(k);
plot(ums/u1,(y-height)/height,'->','Color','#9364CD','Linewidth',2,'MarkerIndices',1:3:length(y1),'MarkerSize',9)

load('shape_triangle9.mat')
y1 = (y-height)/height;
k = y1 >= 0;
k = find(k,1,'first');
u1 = ums(k);
plot(ums/u1,(y-height)/height,'-v','Color','#E32017','Linewidth',2,'MarkerIndices',1:3:length(y1),'MarkerSize',9)

load('shape_semicircle.mat')
y1 = (y-height)/height;
k = y1 >= 0;
k = find(k,1,'first');
u1 = ums(k);
plot(ums/u1,(y-height)/height,'-o','Color','#E32017','Linewidth',2,'MarkerIndices',1:3:length(y1),'MarkerSize',9)

load('shape_trapezium.mat')
y1 = (y-height)/height;
k = y1 >= 0;
k = find(k,1,'first');
u1 = ums(k);
plot(ums/u1,(y-height)/height,'-d','Color','#9364CD','Linewidth',2,'MarkerIndices',1:3:length(y1),'MarkerSize',9)

load('shape_blade.mat')
y1 = (y-height)/height;
k = y1 >= 0;
k = find(k,1,'first');
u1 = ums(k);
plot(ums/u1,(y-height)/height,'-+','Color','#003688','Linewidth',2,'MarkerIndices',1:3:length(y1),'MarkerSize',9)

yline(0,'--','Color','k','Linewidth',2)
set(gcf,'position',[160 280 800 600])
set(gca,'Ylim',[-1 0.2])
set(gca,'Xlim',[0 1.5])
set(gca,'Fontn','Times','FontSize',18,'LineWidth',2)
ylabel('$y^+/k^+$','Interpreter','Latex','FontSize',26)
xlabel('$\bar{u}_g^+/\bar{u}_{g,y=0}^+$','Interpreter','Latex','FontSize',26)
box on