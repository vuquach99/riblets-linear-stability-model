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

figure(1),clf
set(gcf,'position',[1 1 1150 380])
set(gcf,'PaperPositionMode','auto')

subplot(1,2,1),hold on
set(gca,'Xlim',[10 500])
set(gca,'Ylim',[0 .2])
set(gca,'xscale','log')
set(gca,'Fontn','Times','FontSize',20,'linew',2)
set(gca,'YTick',[0:.1:.2])
xlabel('$\lambda_x^+$','interpreter','latex','Pos',[250 -.0107 0]);
ylabel('$\sigma_I^+$','interpreter','latex','Pos',[9 .15 0]);
%xlabel('                \lambda_x^+','Ver','ba','Hor','l','FontAngle','italic');
%ylabel('          \sigma_I^+','Ver','t','Hor','l','FontAngle','italic');
box on

subplot(1,2,2),hold on
set(gca,'Xlim',[10 500])
set(gca,'Ylim',[0 11.5])
set(gca,'xscale','log')
set(gca,'Fontn','Times','FontSize',20,'linew',2)
set(gca,'YTick',[0:5:15])
xlabel('$\lambda_x^+$','interpreter','latex','Pos',[250 -.55 0]);
ylabel('$c_R^+$','interpreter','latex','Pos',[9 7.5 0]);
%xlabel('                \lambda_x^+','Ver','ba','Hor','l','FontAngle','italic');
%ylabel('     c_R^+','Ver','t','Hor','l','FontAngle','italic');
box on

lnsty(1,:)='o ';lncol(1,:)='k';lnwth(1)=2;mrkcol(1,:)=[0 0 0];mrksz(1)= 6;
lnsty(2,:)='- ';lncol(2,:)='k';lnwth(2)=2;mrkcol(2,:)=[1 1 1];mrksz(2)= 1;
%lnsty(3,:)='o ';lncol(3,:)='k';lnwth(3)=2;mrkcol(3,:)=[ 1  1  1];mrksz(3)=11;
%lnsty(4,:)='o ';lncol(4,:)='k';lnwth(4)=.5;mrkcol(4,:)=[.3 .3 .3];mrksz(4)= 6;
lnsty(4,:)='^ ';lncol(4,:)='k';lnwth(4)=2;mrkcol(4,:)=[1 1 1];mrksz(4)=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jRe=[2 4 1]%4%[1:NR]%2%NR:-1:1%
  Rt=Rtt(jRe)
  for jK=1:NK
    Kwpp=Kwp(jK);
    fname=['ribstab_Rt' num2str(Rt) '_Kw' num2str(Kwpp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    maxeigvl=imag(maxeigvl)/utau/Rt;
    eval(['ii=find(lxp>lxpmin'  num2str(Rt) '(jK));'])
    if jRe==3 & Kwpp==4
      ii=[ii(1:end-5) ii(end-3:end)];
    end
    if jRe==1
      subplot(1,2,1)
      if Kwpp==2
        mylogplot(lxp(ii),maxeigvl(ii), ...
             12,.01, ...
             'o',6,'Color','k','MarkerFaceColor','k','Linew',2)
      else
        mylogplot(lxp(ii),maxeigvl(ii), ...
             round(38/log(10/Kwpp^(.2)))+2,.1, ...%15, ...%
             'o',6,'Color','k','MarkerFaceColor','k','Linew',2)
      end
      subplot(1,2,2)
      mylogplot(lxp(ii),Vmaxeig (ii), ...
           round(25/log(10/Kwpp^(.2)))+2,.1, ...
           'o',6,'Color','k','MarkerFaceColor','k','Linew',2)
    elseif jRe==4
      subplot(1,2,1)
      if Kwpp==2
        mylogplot(lxp(ii),maxeigvl(ii), ...
             10,.01, ...
             '^',11,'Color','k','MarkerFaceColor','w','Linew',2)
      else
        mylogplot(lxp(ii),maxeigvl(ii), ...
             round(38/log(10/Kwpp^(.2))),.1, ...%15, ...%
             '^',11,'Color','k','MarkerFaceColor','w','Linew',2)
      end
      subplot(1,2,2)
      mylogplot(lxp(ii),Vmaxeig (ii), ...
           round(25/log(10/Kwpp^(.2))),.1, ...
           '^',11,'Color','k','MarkerFaceColor','w','Linew',2)
    else
%colcode=[(jK-1)/(NK-1) 0 1-(jK-1)/(NK-1)];%
%jRe,lnsty(jRe,:)
      subplot(1,2,1)
      plot(lxp(ii),maxeigvl(ii), ...%plot(lxp,maxeigvl, ...%
           'LineStyle', lnsty(jRe,:), ...
           'Color'    , lncol(jRe,:), ...%'Color'    , colcode, ...%
           'LineWidth', lnwth(jRe  ), ...
           'MarkerSiz', mrksz(jRe  ), ...
           'MarkerFac',mrkcol(jRe,:))
      subplot(1,2,2)
      plot(lxp(ii),Vmaxeig(ii), ...%plot(lxp,Vmaxeig, ...%
           'LineStyle', lnsty(jRe,:), ...
           'Color'    , lncol(jRe,:), ...%'Color'    , colcode, ...%
           'LineWidth', lnwth(jRe  ), ...
           'MarkerSiz', mrksz(jRe  ), ...
           'MarkerFac',mrkcol(jRe,:))
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
annotation(gcf,'arrow',[0.2583 0.2583],[0.3842 0.68],...
           'Color',[.3 .3 .3],'HeadLength',15,'HeadWidth',15,'LineWidth',4);
subplot(1,2,2)
annotation(gcf,'arrow',[0.7826 0.8539],[0.4447 0.6895],...
           'Color',[.3 .3 .3],'HeadLength',15,'HeadWidth',15,'LineWidth',4);
print -depsc2 new_ampl_convvel_exp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kwp=[2 3 4 5 6 8 10 12]';
NK=size(Kwp,1);
lnsty(1,:)='--';lncol(1,:)='k';lnwth(1)=2;mrkcol(1,:)=[1 1 1];mrksz(1)= 1;
lnsty(2,:)='- ';lncol(2,:)='k';lnwth(2)=2;mrkcol(2,:)=[1 1 1];mrksz(2)= 1;
lnsty(3,:)='o ';lncol(3,:)='k';lnwth(3)=2;mrkcol(3,:)=[0 0 0];mrksz(3)= 6;
lnsty(4,:)='^ ';lncol(4,:)='k';lnwth(4)=2;mrkcol(4,:)=[1 1 1];mrksz(4)=13;


figure(2),clf,hold on
set(gcf,'position',[1 1 500 380])
set(gcf,'PaperPositionMode','auto')
%set(gcf,'position',[1 1 500 380])
%set(gcf,'PaperPositionMode','auto')
%set(gca,'position',[.13 .16 .77 .70])
set(gca,'Xlim',[0 12])
set(gca,'Ylim',[0 .2])
set(gca,'Fontn','Times','FontSize',20,'Linew',2)
set(gca,'XTick',[0:6:24])
set(gca,'YTick',[0:.1:.8])
%xlabel('           L_w^+','Ver','ba','Hor','l','FontAngle','italic');
%ylabel('  \sigma_I^+_{ max}','Ver','t','Hor','c','FontAngle','italic');
xlabel('$L_w^+$','interpreter','latex','Pos',[9 -.005 0]);
ylabel('$\sigma_{I\,max}^+$','interpreter','latex','Pos',[-.5 .15 0]);
box on
maxampl=zeros(NK,NR);
for jRe=[1 2 4 3]
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
       'LineWidth', lnwth(jRe  ), ...
       'MarkerSiz', mrksz(jRe  ), ...
       'MarkerFac',mrkcol(jRe,:))
end
print -depsc2 new_maxampl_exp
