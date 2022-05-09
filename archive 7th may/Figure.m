clear
clc
Nosmod=[512];
% Nosmod=256;
nosmod=Nosmod;
% Lwp=[10]';%[2 4 6 8 10 12]';
% Rtt=[180 550 1000]';
Rtt=[550];
% NK=size(Lwp,1);
% NR=size(Rtt,1);
% s=[1 4 6 8 10 10.5 11 11.5 11.75 12 14 16 18 20 25 30 40 50]'; % riblet spacing in wall units [2 4 6 8 10 12 14]';
% % s=[4 6 8 10 12 16 20 30 40 50]';
% Lwp=((pi/256)^(1/3))*s;

lxpmin180 =[14 14 15 15 14 14 14 14 14];
lxpmin400 =[14 14 14 14 14 14 14 14 14];
lxpmin550 =[14 14 14 14 14 14 14 14 14];
lxpmin1000 =[14 14 14 14 14 14 14 14 14];
lxpmin2000=[14 18 18 18 18 18 14 14 14];
lxpmin10000=[14 18 18 18 18 18 14 14 14];
zi=sqrt(-1);

% % Plotting parameters
% lnsty(1,:)='- ';lncol(1,:)='r';lnwth(1)=2;mrkcol(1,:)=[0 0 0];mrksz(1)= 6;
% lnsty(2,:)='--';lncol(2,:)='b';lnwth(2)=4;mrkcol(2,:)=[1 1 1];mrksz(2)= 1;
% lnsty(3,:)='-.';lncol(3,:)='m';lnwth(3)=3;mrkcol(3,:)=[ 1  1  1];mrksz(3)=11;
% lnsty(4,:)='- ';lncol(4,:)='k';lnwth(4)=2;mrkcol(4,:)=[1 1 1];mrksz(4)=10;
% lnsty(5,:)='--';lncol(5,:)='r';lnwth(5)=3;mrkcol(5,:)=[1 1 1];mrksz(4)=6;

% ak2=alp0.^2;
mrk{1}='k-o'; mrk{2}='k-s'; mrk{3}='k-^'; mrk{4}='k-d';
clr(1,:)='k'; clr(2,:)='k'; clr(3,:)='k'; clr(4,:)='k';
line{1} = 'k-'; line{2} = 'k--'; line{3} = 'k:';

s=[10 14 16 18 20 30 40 50]';%[1 2 6 10 11 14 16 50 60 70 80]';%[1 2 6 8 10 11 12 14 16 18 20 25 30 40 50]';
% Lwp=((pi/256)^(1/3))*s;
Lwp=((pi/256)^(1/3))*s;
% flowtype='viscous';
%% Triangular riblets
% F1=0.004448936237952;
% F2=0.020842780639986;
%% Trapezoidal riblets
% F1=0.007178406796948;
% F2=0.027606404326489;
% %% Blade
% F1=0.0083;
% F2=0.0229;
% %% Permeabilities
% Lwp=F1^(1/3)*s;
% Lsp=F2^(1/2)*s;

%% Plot sigma I (amplification) vs wavelength
 figure
% ,clf,hold on
hold on
% subplot(1,2,1)
% figure
set(gcf,'position',[1 1 575 380])
set(gcf,'PaperPositionMode','auto')
        
Ax = gca;
Color = Ax.ColorOrder;
for jRe=1:length(Rtt)
    Rt=Rtt(jRe);
    jnosmod=Nosmod(jRe);
    marker=mrk{jRe};
    lne=line{1};
    color=clr(jRe,:);
    for jK=1:length(Lwp)
        Lwpp=Lwp(jK);
%         fname=['results/nosmod_' num2str(jnosmod) '_Rt_' num2str(Rt) '_Lwp_' num2str(Lwpp) '_' flowtype '.mat'];
        fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat']
        load(fname)
        nx=length(lxp);
%          utau=Rt/Re;
        for p=1:nx;
        imag_eigval(p)=max(imag(eigvals(:,p)))/utau/Rt; %.*Re./Rt.^2;
%         imag_eigval(p)=max(imag(eigvalsPos{p}))/utau/Rt; %.*Re./Rt.^2;
        end
        hold on
        n = length(Lwp);
%         plot(lxp(4:end),imag_eigval(4:end), marker, 'LineWidth', 1, 'Color', [(jK-1)*1/n', 0, 1-(jK-1)*1/n])
        plot(lxp,imag_eigval, lne, 'LineWidth', 2, 'Color', [(jK-1)*1/n', 0, 1-(jK-1)*1/n])
%         plot(lxp,maxeigvl, marker, 'LineWidth', 2, 'Color', [(jK-1)*1/n', 0, 1-(jK-1)*1/n])
        title('Amplification vs Wave-length')
        set(gca,'xscale','log')
        set(gca,'Fontn','Times','FontSize',14,'LineWidth',2)
                set(gca,'YTick',[0:.1:.2])
%         xlim([10 500])
%         ylim([0 0.2])
        set(gca,'Xlim',[10 700])
        set(gca,'Ylim',[0 0.2])
        xlabel('\lambda_x^+','FontAngle','italic'); %'Ver','ba','Hor','l',
        ylabel('\sigma_{Im}^+','FontAngle','italic'); %'Ver','t','Hor','l',
        box on
%         annotation(gcf,'arrow',[0.4583 0.4583],[0.3842 0.68],...
%             'Color',[.3 .3 .3],'HeadLength',15,'HeadWidth',15,'LineWidth',4);
        print -depsc2 ampl_AVERT_1
    end
%     lgn=legend('L_w^+=2','L_w^+=6','L_w^+=10', 'L_w^+=14');
%     title(lgn, 'Wall Normal Permeability Values')
%     set(lgn,'EdgeColor','k','linew',2,'Position',[0.68 0.42 0.25 0.25])
end

%% Plot maximum amplification for each wall normal permeability
figure,clf,hold on
% subplot(1,2,2)
hold on
% figure
Ax = gca;
% Color = Ax.ColorOrder;
% s=[1 6 8 10 10.5 11 12 14 16 18 20 25 30 40 50]'; % riblet spacing in wall units [2 4 6 8 10 12 14]';
% Lwp=((pi/256)^(1/3))*s;

for jRe=1:length(Rtt)
    Rt=Rtt(jRe);
    jnosmod=Nosmod(jRe);
    marker=mrk{3}; %mrk{jRe};
    for kk=1:length(Lwp)
        Lwpp=Lwp(kk);
%         fname=['results/nosmod_' num2str(jnosmod) '_Rt_' num2str(Rt) '_Lwp_' num2str(Lwpp) '_' flowtype '.mat'];
        fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(jnosmod) '.mat']
        load(fname)
        nx=length(lxp);
%         utau=Rt/Re;
        for p=1:nx;
% %         [maxi, zz]=max(imag(eigvals(:,p)));%/utau/Rt;%.*Re./Rt.^2;
%         [real_e, mm]=min(real(eigvals(:,p)));
%         if real_e<0;
%             eigvals(mm,p)=NaN;
%         end
%          [maxi, zzz]=max(imag(eigvals(:,p)))/utau/Rt;%.*Re./Rt.^2;
%         alleigval_imag(p)=maxi/utau/Rt;
        alleigval_imag(p)=max(imag(eigvals(:,p)))/utau/Rt;
%         alleigval_imag(p)=max(imag(eigvalsPos{p}))/utau/Rt;
        end
%         maxampl(kk,jRe)=max(alleigval_imag(4:end));
         maxampl(kk,jRe)=max(alleigval_imag);
    end
    n = length(Lwp);    
    Plot(jRe) = plot(Lwp, maxampl(:,jRe), 'k:' , 'LineWidth', 2, 'MarkerFaceColor', [1 1 1]);
    hold on
    for kk = 1:length(Lwp)
        plot(Lwp(kk),maxampl(kk,jRe), marker, 'MarkerSize',6, 'MarkerFaceColor', [(kk-1)*1/n', 0, 1-(kk-1)*1/n])
    %'LineWidth', 4,'Color',
    end
%  plot([Lwp],maxampl(kk,jRe), marker, 'LineWidth', 1, 'Color',[(kk-1)*1/n', 0, 1-(kk-1)*1/n])
%     set(gcf,'position',[1 600 700 380])
%     set(gcf,'PaperPositionMode','auto')
%     set(gca,'position',[0.10 0.1632 0.56 0.7618])
%     set(gca,'Xlim',[0 0.2])
%     set(gca,'Ylim',[0 .2])
leg{jRe} = ['Re_{\tau}=', num2str(Rt)];
end
% title('Max Amplification vs Wall-normal Permeability')
lgn=legend(Plot, {'Triangular','Trapezoidal','Blade'});
title(lgn,'Riblet Shape')
% lgn=legend(Plot, leg);

    set(gca,'Fontn','Times','FontSize',14,'LineWidth',2)
    set(gca,'XTick',[0:6:24])
    set(gca,'YTick',[0:.2:.8])
    xlabel('L_w^+','FontAngle','italic'); %'Ver','ba','Hor','l',
    ylabel(' \sigma_I^+_{ max}','FontAngle','italic'); %'Ver','t','Hor','c'
    box on

% set(lgn,'EdgeColor','k','linew',2,'Position',[0.68 0.42 0.25 0.25])

print -depsc2 ampl_AVERT_2

% %% Comparison - Plot sigma I (amplification) vs wavelength
% figure
% hold on
% % ,clf,hold on
% % subplot(1,2,1)
% % figure
% set(gcf,'position',[1 1 575 380])
% set(gcf,'PaperPositionMode','auto')
%        
% Ax = gca;
% Color = Ax.ColorOrder;
% % ftypes={ 'inviscid50','pressure50'}
% ftypes={'only_shear'};
% % Rt=Rtt;
% % flowtype='pressure'
% % fname1=['results9thJulyGood/ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat']
% % fname2=['results9thJulyGood/nosmod_' num2str(Nosmod) '_Rt_' num2str(Rtt) '_Lwp_' num2str(Lwpp) '_' flowtype '.mat']   
% % fnames={fname1,fname}
% for c=1:length(ftypes)%c=1:length(fnames) 
%     flowtype=ftypes{c};
%     marker=line{c};
%     color=clr(c,:);
%     Rt=Rtt;
%     nosmod=Nosmod;
% %     fname=fnames{c}
%     for jK=1:length(Lwp)
%         Lwpp=Lwp(jK);
% %         jKnosmod=Nosmod(jK);
%         fname=['results/nosmod_' num2str(Nosmod) '_Rt_' num2str(Rtt) '_Lwp_' num2str(Lwpp) '_' flowtype '.mat'];
% %         fname1=['results9thJulyGood/ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat']
%         load(fname)
%         nx=length(lxp);
%          utau=Rt/Re;
%         for p=1:nx;
% %         imag_eigval(p)=max(imag(eigvals(:,p)))/utau/Rt; %.*Re./Rt.^2;
%         imag_eigval(p)=max(imag(eigvalsPos{p}))/utau/Rt; %.*Re./Rt.^2;
%         end
%         hold on
%         n = length(Lwp);
% %         plot(lxp(4:end),imag_eigval(4:end), marker, 'LineWidth', 1, 'Color', [(jK-1)*1/n', 0, 1-(jK-1)*1/n])
%         plot(lxp,imag_eigval, marker, 'LineWidth', 1, 'Color', [(jK-1)*1/n', 0, 1-(jK-1)*1/n])
%         title('Amplification vs Wave-length')
%         set(gca,'xscale','log')
%         set(gca,'Fontn','Times','FontSize',14,'LineWidth',2)
%         %         set(gca,'YTick',[0:.1:.2])
% %         xlim([10 500])
% %         ylim([0 0.2])
% %         set(gca,'Xlim',[10 700])
% %         set(gca,'Ylim',[0 0.2])
%         xlabel('\lambda_x^+','FontAngle','italic'); %'Ver','ba','Hor','l',
%         ylabel('\sigma_I^+','FontAngle','italic'); %'Ver','t','Hor','l',
%         box on
%         % annotation(gcf,'arrow',[0.4583 0.4583],[0.3842 0.68],...
%         %     'Color',[.3 .3 .3],'HeadLength',15,'HeadWidth',15,'LineWidth',4);
%         print -depsc2 ampl_AVERT_1
%     end
% %     lgn=legend('L_w^+=2','L_w^+=6','L_w^+=10', 'L_w^+=14');
% %     title(lgn, 'Wall Normal Permeability Values')
%     % set(lgn,'EdgeColor','k','linew',2,'Position',[0.68 0.42 0.25 0.25])
% end
% % 
% %% Comparison Plot maximum amplification for each wall normal permeability
% % figure,clf,hold on
% % subplot(1,2,2)
% figure
% hold on
% Ax = gca;
% Color = Ax.ColorOrder;
% ftypes={'pressure50'};
% % s=[1 6 8 10 10.5 11 12 14 16 18 20 25 30 40 50]'; % riblet spacing in wall units [2 4 6 8 10 12 14]';
% % Lwp=((pi/256)^(1/3))*s;
% for d=1:length(ftypes)
%     marker=mrk{d};
%     flowtype=ftypes{d};
% %     for jRe=[1:length(Rtt)]
% %         Rt=Rtt(jRe);
% %         jnosmod=Nosmod(jRe);
%         for kk=1:length(Lwp)
%             Lwpp=Lwp(kk);
%             fname=['results/nosmod_' num2str(Nosmod) '_Rt_' num2str(Rtt) '_Lwp_' num2str(Lwpp) '_' flowtype '.mat'];
%             load(fname)
%             nx=length(lxp);
%             utau=Rt/Re;
%             for p=1:nx;
%     % %         [maxi, zz]=max(imag(eigvals(:,p)));%/utau/Rt;%.*Re./Rt.^2;
%     %         [real_e, mm]=min(real(eigvals(:,p)));
%     %         if real_e<0;
%     %             eigvals(mm,p)=NaN;
%     %         end
%     %          [maxi, zzz]=max(imag(eigvals(:,p)))/utau/Rt;%.*Re./Rt.^2;
%     %         alleigval_imag(p)=maxi/utau/Rt;
%     %         alleigval_imag(p)=max(imag(eigvals(:,p)))/utau/Rt;
%             alleigval_imag(p)=max(imag(eigvalsPos{p}))/utau/Rt;
%             end
%     %         maxampl(kk,jRe)=max(alleigval_imag(4:end));
%              maxampl(kk,d)=max(alleigval_imag);
%         end
%         n = length(Lwp);    
%         Plot(d) = plot(Lwp, maxampl(:,d), marker, 'LineWidth', 2, 'MarkerFaceColor', [1 1 1]);
%         for kk = 1:length(Lwp)
%             plot(Lwp(kk),maxampl(kk,d), marker, 'MarkerSize',6, 'MarkerFaceColor', [(kk-1)*1/n', 0, 1-(kk-1)*1/n])
%         %'LineWidth', 4,'Color',
%         end
%     %  plot([Lwp],maxampl(kk,jRe), marker, 'LineWidth', 1, 'Color',[(kk-1)*1/n', 0, 1-(kk-1)*1/n])
%     %     set(gcf,'position',[1 600 700 380])
%     %     set(gcf,'PaperPositionMode','auto')
%     %     set(gca,'position',[0.10 0.1632 0.56 0.7618])
%     %     set(gca,'Xlim',[0 0.2])
%     %     set(gca,'Ylim',[0 .2])
%     leg{d} = ftypes{d};
% %     end
% end
% title('Max Amplification vs Wall-normal Permeability')
% % lgn=legend(Plot, {'Re_{\tau}=180','Re_{\tau}=550','Re_{\tau}=1000','Re_{\tau}=2000'});
% lgn=legend(Plot, leg);
% 
%     set(gca,'Fontn','Times','FontSize',14,'LineWidth',2)
% %     set(gca,'XTick',[0:6:24])
% %     set(gca,'YTick',[0:.2:.8])
%     xlabel('L_w^+','FontAngle','italic'); %'Ver','ba','Hor','l',
%     ylabel(' \sigma_I^+_{ max}','FontAngle','italic'); %'Ver','t','Hor','c'
%     box on
% 
% set(lgn,'EdgeColor','k','linew',2,'Position',[0.68 0.42 0.25 0.25])
% 
% print -depsc2 ampl_AVERT_2

%% Plot Orr-Sommerfeld spectrum (cI vs cR)
% % figure,clf,hold on
% % maxampl=zeros(NK,NR);
% % omega_imag=zeros(nosmod+1,nx);
% % omega_real=zeros(nosmod+1,nx);
% % for jRe=[1:NR]
% %     Rt=Rtt(jRe)
% %     for z=[1:NK]
% %         Lwpp=Lwp(z);
% % %       subplot(2,3,z)
% %         fname=['results/turb_ribstab_Rt' num2str(Rt) '_Lw_' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat'];
% %         load(fname)
% %         omega_imag=imag(eigvals);%/utau/Rt;%wall units
% %         omega_real=real(eigvals);%/utau/Rt;%wall units
% %         
% %         hold on
% %         plot(omega_real, omega_imag, 'o', 'LineWidth', 2) %, 'Color', lncol(jK,:),
% %   %     plot(c_real./alp0, c_imag./alp0, 'o', 'LineWidth', 2) %, 'Color', lncol(jK,:),
% %         %         'LineStyle', lnsty(jRe,:), ...
% %         %         'Color'    , lncol(jRe,:), ...
% %         %         'LineWidth', lnwth(jRe  ))
% %         % set(gca,'Xlim',[10 500])
% %         % set(gca,'Ylim',[0 .2]) ./alp0
% %         title('O-S Spectrum Viscous du/dy=0, Re=185 ')
% %         set(gca,'Fontn','Times','FontSize',18,'LineWidth',2)
% %         %         set(gca,'YTick',[0:.1:.2])
% %         xlabel('w_R^+','FontAngle','italic'); %'Ver','ba','Hor','l',
% %         ylabel('w_I^+','FontAngle','italic'); %'Ver','t','Hor','l',
% %         box on
% % %         legenda=['L_w^+=2','L_w^+=6','L_w^+=10', 'L_w^+=14']
% % %         lgn=legenda(z);
% % %         title(lgn, 'Wall Normal Permeability Values');
% %     end
% %     
% % end