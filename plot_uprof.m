clear
close all

% inflectional velocity profile
x = linspace(0,0.9);
y1 = (1/7)*tan(3*x-1.2925);

figure(1)
plot(x,y1,'LineWidth',2,'Color','b')
yline(0,'--','LineWidth',2,'Color','k')
xline(0,'-','LineWidth',2,'Color','k')
yline(-0.5,'-','LineWidth',2,'Color','k')
set(gcf,'position',[160 280 800 600])
axis equal
set(gca,'Xlim',[-0.1 0.87])
set(gca,'Ylim',[-0.6 0.5])
set(gca,'xtick',[])
set(gca,'ytick',[])
box on

y2 = exp(4*x-3)-0.55;
figure(2)
plot(x,y2,'LineWidth',2,'Color','b')
yline(0,'--','LineWidth',2,'Color','k')
xline(0,'-','LineWidth',2,'Color','k')
yline(-0.5,'-','LineWidth',2,'Color','k')
set(gcf,'position',[160 280 800 600])
axis equal
set(gca,'Xlim',[-0.1 0.87])
set(gca,'Ylim',[-0.6 0.5])
set(gca,'xtick',[])
set(gca,'ytick',[])
box on