nosmod=1536;
Kwp=[2 4 6 8 10 12 14]';
Rtt=[185 550 950 2000]';
%Lwp=[1.5 3 6 12]';
%Rtt=[185 550 950]';
NK=size(Kwp,1);
NR=size(Rtt,1);

eval(['load basefuns2_' num2str(nosmod) ';'])
nm=round(N/2);

plist1=[1 35 50 60 66:4:86 94:8:110 117:6:nm];
plist2=[1 35 50 60 70:8:86 94:8:110 117:6:nm];

lnsty(1,:)='--';lncol(1,:)='k';lnwth(1)=2;
lnsty(2,:)='-.';lncol(2,:)='k';lnwth(2)=2;
lnsty(3,:)='- ';lncol(3,:)='k';lnwth(3)=2;
lnsty(4,:)='--';lncol(4,:)='k';lnwth(4)=2;
lnsty(5,:)='- ';lncol(5,:)='k';lnwth(5)=1;
lnsty(6,:)='- ';lncol(6,:)='k';lnwth(6)=1;
lnsty(7,:)='no';lncol(7,:)='k';lnwth(7)=1;
lnsty(8,:)='- ';lncol(8,:)='k';lnwth(8)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1),clf
set(gcf,'position',[1 400 1150 380])
set(gcf,'PaperPositionMode','auto')

subplot(1,2,1),hold on
set(gca,'Xlim',[0 20])
set(gca,'Ylim',[0 1])
set(gca,'Fontn','Times','FontSize',20,'linew',2)
set(gca,'XTick',[0:10:80])
set(gca,'YTick',[0: 1: 5])
xlabel('$y^+$','interpreter','latex','Pos',[15 -.03 0]);
ylabel('$\left| \widehat{u} \right|$','interpreter','latex','Pos',[-.5 .5 0]);
%xlabel('               y^+','Ver','ba','Hor','l','FontAngle','italic');
%%ylabel('|u|','Ver','ba','Hor','c','FontAngle','italic');
%ylabel('AAA','Ver','ba','Hor','c','FontAngle','italic');
box on

subplot(1,2,2),hold on
set(gca,'Xlim',[0 20])
set(gca,'Ylim',[0 1])
set(gca,'Fontn','Times','FontSize',20,'linew',2)
set(gca,'XTick',[0:10:80])
set(gca,'YTick',[0: 1: 5])
xlabel('$y^+$','interpreter','latex','Pos',[15 -.03 0]);
ylabel('$\left| \widehat{v} \right|$','interpreter','latex','Pos',[-.5 .5 0]);
%xlabel('               y^+','Ver','ba','Hor','l','FontAngle','italic');
%%ylabel('|v|','Ver','ba','Hor','c','FontAngle','italic');
%ylabel('BBB','Ver','ba','Hor','c','FontAngle','italic');
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jRe=2
Rt=Rtt(jRe)
for jK=[2 3 4 7]
  Kwpp=Kwp(jK);
  fname=['ribstab_Rt' num2str(Rt) '_Kw' num2str(Kwpp) '_Ny' num2str(nosmod) '.mat'];
  load(fname)
  [Amax,i]=max(imag(maxeigvl));
%  lxpi=lxp(i)
  alp=2*pi*Rt/lxp(i);
  uuu=abs(D1*maxeigvc(:,i)/alp);
  vvv=abs(D0*maxeigvc(:,i));
  vvv=.5*(vvv(end:-1:nm+2)/vvv(end)+vvv(1:nm)/vvv(1));
  uuu=.5*(uuu(end:-1:nm+2)/uuu(end)+uuu(1:nm)/uuu(1));
  subplot(1,2,1)
  plot(Rt*(1-y(1:nm)),uuu, ...
       'LineStyle', lnsty(jK,:), ...
       'Color'    , lncol(jK,:), ...
       'LineWidth', lnwth(jK  ))
  subplot(1,2,2)
  plot(Rt*(1-y(1:nm)),vvv, ...
       'LineStyle', lnsty(jK,:), ...
       'Color'    , lncol(jK,:), ...
       'LineWidth', lnwth(jK  ))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(1,2,1)
  plot(Rt*(1-y(plist1)),uuu(plist1),'.', ...
       'Color'    , lncol(jK,:), ...
       'LineWidth', lnwth(jK  ), ...
       'MarkerSiz',24, ...
       'MarkerFac','k')
  subplot(1,2,2)
  plot(Rt*(1-y(plist2)),vvv(plist2),'.', ...
       'Color'    , lncol(jK,:), ...
       'LineWidth', lnwth(jK  ), ...
       'MarkerSiz',24, ...
       'MarkerFac','k')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%subplot(1,2,1)
%annotation(gcf,'arrow',[0.4513 0.4487],[0.2947 0.1316],...
%           'Color',[.3 .3 .3],'HeadLength',15,'HeadWidth',15,'LineWidth',4);
%subplot(1,2,2)
%annotation(gcf,'arrow',[0.7678 0.6661],[0.7816 0.1289],...
%           'Color',[.3 .3 .3],'HeadLength',15,'HeadWidth',15,'LineWidth',4);

print -depsc2 new_large_mode
