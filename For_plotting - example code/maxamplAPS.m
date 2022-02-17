nosmod=1536%256%
Kwp=[2 4 6 8 10 12]';
Rtt=[185 550 950 2000]';
NK=size(Kwp,1);
NR=size(Rtt,1);

lxpmin185 =[14 14 15 15 14 14 14 14 14];
lxpmin550 =[14 14 14 14 14 14 14 14 14];
lxpmin950 =[14 14 14 14 14 14 14 14 14];
lxpmin2000=[14 18 18 18 18 18 14 14 14];

zi=sqrt(-1);

lnsty(1,:)='- ';lncol(1,:)='r';lnwth(1)=2;mrkcol(1,:)=[0 0 0];mrksz(1)= 6;
lnsty(2,:)='--';lncol(2,:)='b';lnwth(2)=4;mrkcol(2,:)=[1 1 1];mrksz(2)= 1;
lnsty(3,:)='-.';lncol(3,:)='m';lnwth(3)=3;mrkcol(3,:)=[1 1 1];mrksz(3)=11;
lnsty(4,:)='- ';lncol(4,:)='k';lnwth(4)=2;mrkcol(4,:)=[1 1 1];mrksz(4)=10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kwp=[2 3 4 5 6 8 10 12]';
NK=size(Kwp,1);


figure(2),clf,hold on
%set(gcf,'position',[1 1 460 360])
%set(gcf,'PaperPositionMode','auto')
%set(gca,'position',[0.13 0.15 0.775 0.755])
set(gca,'Xlim',[0 12])
set(gca,'Ylim',[0 .2])
set(gca,'Fontn','Times','FontSize',26,'LineWidth',2)
set(gca,'XTick',[0:6:24])
set(gca,'YTick',[0:.1:.8])
xlabel('$L_w^+$','interpreter','latex','Position',[9 -0.007 0]);
ylabel('$\sigma^+_{I\,max}$','interpreter','latex','Position',[-.5 0.15 0]);
box on
maxampl=zeros(NK,NR);
for jRe=[1 2 4]
  Rt=Rtt(jRe)
  for jK=1:NK
    Kwpp=Kwp(jK);
    fname=['ribstab_Rt' num2str(Rt) '_Kw' num2str(Kwpp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    maxeigvl=imag(maxeigvl)/utau/Rt;
    maxampl(jK,jRe)=max(maxeigvl);
  end
  plot([0;Kwp],[0;maxampl(:,jRe)], ...
       'LineStyle', lnsty(jRe,:), ...
       'Color'    , lncol(jRe,:), ...
       'LineWidth', lnwth(jRe  ))
end

lgn=legend('Re$_{\tau}=185$','Re$_{\tau}=550$','Re$_{\tau}=2000$');
set(lgn, ...
    'interpreter','latex', ...
    'FontSize',24, ...
    'EdgeColor','k', ...
    'linew',2, ...
    'Position',[0.58 0.31 0.25 0.25])
%    'Position',[0.60 0.16 0.25 0.25])

print -depsc2 maxamplAPS

set(gca,'pos',[2 2 1 1])
set(lgn,'Position',[0.68 0.75 0.25 0.25])
print -depsc2 legendAPS
