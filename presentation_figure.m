clear
close all

% Nosmod=[512];
Nosmod=256;
nosmod=Nosmod;
% Lwp=[10]';%[2 4 6 8 10 12]';
% Rtt=[180 550 1000]';
Rtt=[550];
% NK=size(Lwp,1);
% NR=size(Rtt,1);
% s=[1 4 6 8 10 10.5 11 11.5 11.75 12 14 16 18 20 25 30 40 50]'; % riblet spacing in wall units [2 4 6 8 10 12 14]';
% % s=[4 6 8 10 12 16 20 30 40 50]';
% Lwp=((pi/256)^(1/3))*s;

% ak2=alp0.^2;
mrk{1}='k-o'; mrk{2}='k-s'; mrk{3}='k-^'; mrk{4}='k-d';
clr(1,:)='k'; clr(2,:)='k'; clr(3,:)='k'; clr(4,:)='k';
line{1} = 'k-'; line{2} = 'k--'; line{3} = 'k:';

s=[20]';%[1 2 6 10 11 14 16 50 60 70 80]';%[1 2 6 8 10 11 12 14 16 18 20 25 30 40 50]';
% Lwp=((pi/256)^(1/3))*s;
Lwp=((pi/256)^(1/3))*s;
% flowtype='viscous';

% Plot sigma I (amplification) vs wavelength
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
        %fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat']
        fname = ['ribstab_Rt550_Lw2.843_Ny256.mat'] % set inside grooves
        load(fname)
        nx=length(lxp);
%          utau=Rt/Re;
        for p=1:nx;
        imag_eigval(p)=max(imag(eigvals(:,p)))/utau/Rt; %.*Re./Rt.^2;
%         imag_eigval(p)=max(imag(eigvalsPos{p}))/utau/Rt; %.*Re./Rt.^2;
        end
        hold on
        n = length(Lwp);
        plot(lxp,imag_eigval, 'LineWidth', 2)
        fname = ['ribstab_Rt550_Lw4.6134_Ny256.mat'] % set at tips
        load(fname)
        for p=1:nx;
        imag_eigval(p)=max(imag(eigvals(:,p)))/utau/Rt; %.*Re./Rt.^2;
        end
        n = length(Lwp);
        plot(lxp,imag_eigval, 'LineWidth', 2)
        yline(0,'--','LineWidth',2)
        legend('Set inside grooves', 'Set at tips','', 'location', 'Southeast')
        set(gcf,'position',[160 280 800 600])
        title('Amplification vs Wavelength for s = 20, semi-circle')
        set(gca,'xscale','log')
        set(gca,'Fontn','Times','FontSize',18,'LineWidth',2)
        %set(gca,'YTick',[0:.1:.2])
        set(gca,'Xlim',[5 1000])
        %set(gca,'Ylim',[0 0.1])
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
