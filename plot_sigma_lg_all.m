% Plotting script for maximum amplification growth rate vs lg+
% all riblets, at virtual origin
clear
close all

nosmod = 256;
Rt = 550;

%% Calculate groove cross-section for each riblet types
% lg = s*sqrt(pi/8); % semi-circle
% lg = s*sqrt(2/5); % blade
% lg = s*sqrt(0.5 + sqrt(3)/4); % triangle3
% lg = s*sqrt(sqrt(3))/2; % triangle6
% lg = s/2; % triangle9
% lg = s*sqrt(sqrt(3))/2; % trapezium

%% Plot
figure(2)
hold on
% Semi-circle
s1 = [5 7.5 10 12.5 15 17.5 20 30 40 50 60 70 100]';
shape = 'circle1_s';
for jK = 1:size(s1,1)
    sp = s1(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_1(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Blade
s2 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'blade1_s';
for jK = 1:size(s2,1) 
    sp = s2(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_2(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Trangle3
s3 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'triangle31_s';
for jK = 1:size(s3,1) 
    sp = s3(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_3(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Triangle6
s4 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'triangle61_s';
for jK = 1:size(s4,1) 
    sp = s4(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_4(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Triangle9
s5 = [5 7.5 10 12.5 15 17.5 20 30 40 50 60 70 100]';
shape = 'triangle91_s';
for jK = 1:size(s5,1) 
    sp = s5(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_5(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Trapezium
s6 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'trapezium1_s';
for jK = 1:size(s6,1) 
    sp = s6(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_6(1,jK) = imag(Most_unstab)/ut/Rt;
end

plot(s1*sqrt(pi/8),points_1,'-o','LineWidth',2,'Color','b');
plot(s2*sqrt(2/5),points_2,'-o','LineWidth',2,'Color','m');
plot(s3*sqrt(0.5 + sqrt(3)/4),points_3,'-o','LineWidth',2,'Color','k');
plot(s4*sqrt(sqrt(3))/2,points_4,'-o','LineWidth',2,'Color','g');
plot(s5/2,points_5,'-o','LineWidth',2,'Color','r');
plot(s6*sqrt(sqrt(3))/2,points_6,'-o','LineWidth',2,'Color','c');

set(gcf,'position',[160 280 800 600])
set(gca,'Xlim',[0 35])
set(gca,'Fontn','Times','FontSize',18,'LineWidth',2)
xlabel('l_g^+','FontAngle','italic')
ylabel('\sigma_{Im}^+','FontAngle','italic')
legend('semi-circle','blade','triangle3','triangle6','triangle9','trapezium',...
    'location','Southeast','FontSize',18)
title('Maximum amplification vs l_g^+, set at virtual origin')
box on