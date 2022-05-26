% Compare with inviscid code
clear
close all

nosmod = 512;
Rt = 550;
s = [5 10 15 20 30 40 50]';
% shape_0 = 'blade'; % name to appear on plots

%% Amplification growth rate vs wavelength+
figure(1)
hold on
% Inviscid
for jK = 1:size(s,1)
    sp = s(jK);
    Lwp=((pi/256)^(1/3))*sp;
    fname = ['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    nx = length(lxp);
    for p = 1:nx
        imag_eigval(p) = max(imag(eigvals(:,p)))/utau/Rt;
    end
%     most_unstab = imag(Most_unstab)/ut/Rt;
    h(1) = plot(lxp,imag_eigval,'--','LineWidth',2,'Color',...
        [(jK-1)*1/(length(s))',0,1-(jK-1)*1/(length(s))]);
%     plot(Most_lxp,most_unstab,'o','LineWidth', 2,'Color',...
%         [(jK-1)*1/(length(s))',0,1-(jK-1)*1/(length(s))],'MarkerSize', 10)
end
% Set inside grooves
shape = 'circle_s';
for jK = 1:size(s,1)
    sp = s(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    nx = length(lxp);
    for p = 1:nx
        imag_eigval(p) = imag(Max_unstab(p))/ut/Rt;
    end
    most_unstab = imag(Most_unstab)/ut/Rt;
    h(2) = plot(lxp,imag_eigval,'-','LineWidth',2,'Color',...
        [(jK-1)*1/(length(s))',0,1-(jK-1)*1/(length(s))]);
%     plot(Most_lxp,most_unstab,'x','LineWidth', 2,'Color',...
%         [(jK-1)*1/(length(s))',0,1-(jK-1)*1/(length(s))],'MarkerSize', 10)
end
set(gca,'xscale','log')
yline(0,'--','LineWidth',2)
set(gcf,'position',[160 280 800 600])
set(gca,'Xlim',[10 6000])
set(gca,'Ylim',[0 0.2])
set(gca,'Fontn','Times','FontSize',22,'LineWidth',2)
xlabel('$\lambda_x^+$','Interpreter','latex','FontSize',32)
ylabel('$\sigma_I^+$','Interpreter','latex','FontSize',32)
% legend(h([1,2]),{'Set at tips','Set inside grooves'},'location','Southeast','FontSize',18)
% title(sprintf('Amplification vs Wavelength, %s', shape_0))
box on