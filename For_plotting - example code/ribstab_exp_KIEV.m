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

figure(1),clf,hold on
set(gcf,'position',[1 1 460 720])
set(gcf,'PaperPositionMode','auto')
%set(gcf,'position',[1 1 500 800])
%set(gcf,'PaperPositionMode','auto')

lnsty(1,:)='- ';lncol(1,:)='r';lnwth(1)=2;mrkcol(1,:)=[0 0 0];mrksz(1)= 6;
lnsty(2,:)='--';lncol(2,:)='b';lnwth(2)=4;mrkcol(2,:)=[1 1 1];mrksz(2)= 1;
lnsty(3,:)='-.';lncol(3,:)='m';lnwth(3)=3;mrkcol(3,:)=[1 1 1];mrksz(3)=11;
lnsty(4,:)='- ';lncol(4,:)='k';lnwth(4)=2;mrkcol(4,:)=[1 1 1];mrksz(4)=10;

subplot(2,1,1),hold on
set(gca,'Xlim',[10 500])
set(gca,'Ylim',[0 .2])
set(gca,'xscale','log')
set(gca,'Fontn','Times','FontSize',22,'LineWidth',2)
set(gca,'YTick',[0:.1:.2])
xlabel('                \lambda_x^+','Ver','m','Hor','l','FontAngle','italic');
ylabel('          \sigma_I^+','Ver','t','Hor','l','FontAngle','italic');
box on
annotation(gcf,'arrow',[0.45 0.45],[0.67 0.85],...
           'Color','m','HeadLength',15,'HeadWidth',15,'LineWidth',4);

subplot(2,1,2),hold on
set(gca,'Xlim',[10 500])
set(gca,'Ylim',[0 11.5])
set(gca,'xscale','log')
set(gca,'Fontn','Times','FontSize',22,'linew',2)
set(gca,'YTick',[0:5:15])
xlabel('                \lambda_x^+','Ver','m','Hor','l','FontAngle','italic');
ylabel('     c_R^+','Ver','t','Hor','l','FontAngle','italic');
annotation(gcf,'arrow',[0.60 0.85],[0.18 0.35],...
           'Color','m','HeadLength',15,'HeadWidth',15,'LineWidth',4);
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jRe=[1 2 4]
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
    subplot(2,1,1)
    plot(lxp(ii),maxeigvl(ii), ...%plot(lxp,maxeigvl, ...%
         'LineStyle', lnsty(jRe,:), ...
         'Color'    , lncol(jRe,:), ...%'Color'    , colcode, ...%
         'LineWidth', lnwth(jRe  ))
    if jRe==4 & Kwpp==2
      ii=ii(1:2:end);
    end
    subplot(2,1,2)
    plot(lxp(ii),Vmaxeig(ii), ...%plot(lxp,Vmaxeig, ...%
         'LineStyle', lnsty(jRe,:), ...
         'Color'    , lncol(jRe,:), ...%'Color'    , colcode, ...%
         'LineWidth', lnwth(jRe  ))
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1)
set(gca,'position',[0.13 0.58 0.775 0.377])
subplot(2,1,2)
set(gca,'position',[0.13 0.10 0.775 0.377])

print -depsc2 ribstab_exp_KIEV


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kwp=[2 3 4 5 6 8 10 12]';
NK=size(Kwp,1);


figure(2),clf,hold on
set(gcf,'position',[1 1 460 360])
set(gcf,'PaperPositionMode','auto')
set(gca,'position',[0.13 0.15 0.775 0.755])
set(gca,'Xlim',[0 12])
set(gca,'Ylim',[0 .2])
set(gca,'Fontn','Times','FontSize',22,'LineWidth',2)
set(gca,'XTick',[0:6:24])
set(gca,'YTick',[0:.1:.8])
xlabel('           L_w^+','Ver','ba','Hor','l','FontAngle','italic');
ylabel('      \sigma_I^+_{ max}','Ver','t','Hor','l','FontAngle','italic');
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

lgn=legend('Re_{\tau}=185','Re_{\tau}=550','Re_{\tau}=2000');
set(lgn,'EdgeColor','k','linew',2,'Position',[0.58 0.31 0.25 0.25])

print -depsc2 maxampl_KIEV
