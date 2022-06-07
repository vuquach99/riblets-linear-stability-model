% Plotting script for maximum amplification growth rate vs lg+
% at different nosmod's
% Comparing boundary set at tips vs set at v.o.
clear
close all

Rt = 550;
s = [5 7.5 10 12.5 15 17.5 20 30 40 50 60 70 100]';
%% Calculate groove cross-section for each riblet types
% lg = sqrt(pi/8); % semi-circle
% lg = sqrt(2/5); % blade
% lg = sqrt(0.5 + sqrt(3)/4); % triangle3
% lg = sqrt(sqrt(3))/2; % triangle6
% lg = 0.5; % triangle9
% lg = sqrt(sqrt(3))/2; % trapezium

%% Plot
figure(2)
hold on
nosmod = 128;
shape = 'circle1_s';
for jK = 1:size(s,1) 
    sp = s(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_1(1,jK) = imag(Most_unstab)/ut/Rt;
end

nosmod = 256;
shape = 'circle1_s';
for jK = 1:size(s,1)
    sp = s(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_2(1,jK) = imag(Most_unstab)/ut/Rt;
end

nosmod = 512;
shape = 'circle1_s';
for jK = 1:size(s,1)
    sp = s(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_3(1,jK) = imag(Most_unstab)/ut/Rt;
end

plot(s*sqrt(pi/8),points_1,'-x','LineWidth',2,'Color','#00782A');
plot(s*sqrt(pi/8),points_2,'-x','LineWidth',2,'Color','#EE7C0E');
plot(s*sqrt(pi/8),points_3,'-x','LineWidth',2,'Color','#9B0056');

set(gcf,'position',[160 280 800 600])
set(gca,'Xlim',[0 65])
set(gca,'Ylim',[0 0.3])
set(gca,'Fontn','Times','FontSize',30,'LineWidth',2)
xlabel('$\ell_g^+$','Interpreter','latex','FontSize',38)
ylabel('$\sigma_{I}^+$','Interpreter','latex','FontSize',38)
legend('$N = 128$','$N = 256$','$N = 512$',...
    'location','Southeast','Interpreter','latex','FontSize',26)
box on