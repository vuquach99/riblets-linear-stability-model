% Plotting script for maximum amplification growth rate vs lg+
% Comparing boundary set at tips vs set at v.o.
clear
close all

nosmod = 256;
Rt = 550;

%% Calculate groove cross-section for each riblet types
% lg = sqrt(pi/8); % semi-circle
% lg = sqrt(2/5); % blade
lg = sqrt(0.5 + sqrt(3)/4); % triangle3
% lg = sqrt(sqrt(3))/2; % triangle6
% lg = 0.5; % triangle9
% lg = sqrt(sqrt(3))/2; % trapezium

%% Plot
figure(2)
hold on
% Set at tips
s1 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'triangle3_s';
for jK = 1:size(s1,1)
    sp = s1(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_1(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Set at v.o.
s2 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'triangle31_s';
for jK = 1:size(s2,1) 
    sp = s2(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_2(1,jK) = imag(Most_unstab)/ut/Rt;
end

plot(s1*lg,points_1,'--x','LineWidth',2,'Color','k');
plot(s2*lg,points_2,'-o','LineWidth',2,'Color','k');

set(gcf,'position',[160 280 800 600])
set(gca,'Xlim',[0 35])
set(gca,'Ylim',[0 0.15])
set(gca,'Fontn','Times','FontSize',18,'LineWidth',2)
xlabel('l_g^+','FontAngle','italic')
ylabel('\sigma_{Im}^+','FontAngle','italic')
legend('Set at tips','Set at virtual origin','location','Southeast','FontSize',18)
title('Maximum amplification vs l_g^+, triangle3')
box on