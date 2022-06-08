% Plots Orr-Sommerfeld matrix
clear
close all

%% Geometry parameters
s = [5 10 15 20 30 40 50]'; % riblet spacing s+

Rt = 550; % friction Reynolds number
nosmod = 256; % number of modes

%% Main loop, plots Orr-Sommerfeld spectrum
omega_imag = zeros(nosmod+1,size(s,1));
omega_real = zeros(nosmod+1,size(s,1));
most_imag = [];
most_real = [];
figure
hold on

%% First shape
shape = 'circle';
fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(40) '_Ny' num2str(nosmod) '.mat'];
load(fname)
index = find(lxp==Most_lxp);
for jK = 1:size(s,1)
    sp = s(jK)
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    eigvals = eigvals(:,index);
    omega_imag(:,jK) = imag(eigvals)/ut/Rt; % /ut/Rt for channel units
    omega_real(:,jK) = real(eigvals)/ut/Rt;
    plot(omega_real(:,jK),omega_imag(:,jK),'.',...
        'Color',[(jK-1)*1/(length(s))',0,1-(jK-1)*1/(length(s))],'MarkerSize',14)
end

%% Second shape
% shape = 'circle_s';
% fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(40) '_Ny' num2str(nosmod) '.mat'];
% load(fname)
% index = find(lxp==Most_lxp);
% for jK = 1:size(s,1)
%     sp = s(jK)
%     fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
%     load(fname)
%     eigvals = eigvals(:,index);
%     omega_imag(:,jK) = imag(eigvals)/ut/Rt; % /ut/Rt for channel units
%     omega_real(:,jK) = real(eigvals)/ut/Rt;
%     plot(omega_real(:,jK),omega_imag(:,jK),'+',...
%         'Color',[(jK-1)*1/(length(s))',0,1-(jK-1)*1/(length(s))],'MarkerSize',14)
% end

%% Parameters
set(gcf,'position',[160 280 800 600])
set(gca,'Xlim',[-0.5 2.5])
set(gca,'Ylim',[-2 0.5])
% set(gca,'Xlim',[0.3 1])
% set(gca,'Ylim',[-0.3 0.2])
set(gca,'Fontn','Times','FontSize',32,'LineWidth',2)
xlabel('$\omega_r$','Interpreter','latex','FontSize',40)
ylabel('$\omega_i$','Interpreter','latex','FontSize',40)
% legend('','Most amplified mode','FontSize',18)
% title('Orr-Sommerfeld spectrum')
box on