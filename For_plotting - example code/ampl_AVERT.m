clc
clear

nosmod=512%256%
% s=[5 15 25 35 42 50]'
s=[5 10 20 30 40 50]';
Lwp=((pi/256)^(1/3))*s; 
Rtt=[550]';
NK=size(Lwp,1);
NR=size(Rtt,1);

lxpmin185 =[14 14 15 15 14 14 14 14 14];
lxpmin550 =[14 14 14 14 14 14 14 14 14];
lxpmin950 =[14 14 14 14 14 14 14 14 14];
lxpmin2000=[14 18 18 18 18 18 14 14 14];

zi=sqrt(-1);

figure,clf,hold on
%set(gcf,'position',[1 1 575 380])
%set(gcf,'PaperPositionMode','auto')

lnsty(1,:)='- ';lncol(1,:)='r';lnwth(1)=2;mrkcol(1,:)=[0 0 0];mrksz(1)= 6;
lnsty(2,:)='--';lncol(2,:)='b';lnwth(2)=4;mrkcol(2,:)=[1 1 1];mrksz(2)= 1;
lnsty(3,:)='-.';lncol(3,:)='m';lnwth(3)=3;mrkcol(3,:)=[ 1  1  1];mrksz(3)=11;
lnsty(4,:)='- ';lncol(4,:)='k';lnwth(4)=2;mrkcol(4,:)=[1 1 1];mrksz(4)=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jRe=[1]
  Rt=Rtt(jRe)
  for jK=1:NK
    Kwpp=Lwp(jK);
    fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Kwpp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
%     maxeigvl=imag(maxeigvl)/utau/Rt;
    maxeigvl=max(imag(eigvals))/utau/Rt;
%     maxreal=real(maxeigvl)/utau/Rt;
    eval(['ii=find(lxp>lxpmin'  num2str(Rt) '(jK));'])
    if jRe==3 & Kwpp==4
      ii=[ii(1:end-5) ii(end-3:end)];
    end
    plot(lxp(ii),maxeigvl(ii), ...%plot(lxp,maxeigvl, ...%
         'LineStyle', lnsty(jRe,:), ...
         'Color'    , lncol(jRe,:), ...%'Color'    , colcode, ...%
         'LineWidth', lnwth(jRe  ))
%     plot(real(eigvals),imag(eigvals),'.')
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'Xlim',[10 500])
set(gca,'Ylim',[0 .2])
set(gca,'xscale','log')
set(gca,'Fontn','Times','FontSize',20,'LineWidth',2)
set(gca,'YTick',[0:.1:.2])
xlabel('                \lambda_x^+','Ver','ba','Hor','l','FontAngle','italic');
ylabel('          \sigma_I^+','Ver','t','Hor','l','FontAngle','italic');
box on
annotation(gcf,'arrow',[0.4583 0.4583],[0.3842 0.68],...
           'Color',[.3 .3 .3],'HeadLength',15,'HeadWidth',15,'LineWidth',4);

print -depsc2 ampl_AVERT_1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kwp=[2 3 4 5 6 8 10 12]';
% NK=size(Kwp,1);


figure,clf,hold on
%set(gcf,'position',[1 600 700 380])
%set(gcf,'PaperPositionMode','auto')
%set(gca,'position',[0.10 0.1632 0.56 0.7618])
set(gca,'Xlim',[0 12])
set(gca,'Ylim',[0 .2])
set(gca,'Fontn','Times','FontSize',24,'LineWidth',2)
set(gca,'XTick',[0:6:24])
set(gca,'YTick',[0:.2:.8])
xlabel('           L_w^+','Ver','ba','Hor','l','FontAngle','italic');
ylabel('  \sigma_I^+_{ max}','Ver','t','Hor','c','FontAngle','italic');
box on
maxampl=zeros(NK,NR);
for jRe=[1]
  Rt=Rtt(jRe)
  for jK=1:NK
    Kwpp=Lwp(jK);
    fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Kwpp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
%     maxeigvl=imag(maxeigvl)/utau/Rt;
    maxeigvl=max(imag(eigvals))/utau/Rt;
    maxampl(jK,jRe)=max(maxeigvl);
  end
  plot([0;Lwp],[0;maxampl(:,jRe)], ...
       'LineStyle', lnsty(jRe,:), ...
       'Color'    , lncol(jRe,:), ...
       'LineWidth', lnwth(jRe  ))
end

lgn=legend('Re_{\tau}=185','Re_{\tau}=550','Re_{\tau}=2000');
set(lgn,'EdgeColor','k','linew',2,'Position',[0.68 0.42 0.25 0.25])

print -depsc2 ampl_AVERT_2
