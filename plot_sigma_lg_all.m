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

%% At the virtual origin
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

%% At the plane of the tips
% Semi-circle
s7 = [5 7.5 10 12.5 15 17.5 20 30 40 50 60 70 100]';
shape = 'circle_s';
for jK = 1:size(s7,1)
    sp = s7(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_7(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Blade
s8 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'blade_s';
for jK = 1:size(s8,1) 
    sp = s8(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_8(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Trangle3
s9 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'triangle3_s';
for jK = 1:size(s9,1) 
    sp = s9(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_9(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Triangle6
s10 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'triangle6_s';
for jK = 1:size(s10,1) 
    sp = s9(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_10(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Triangle9
s11 = [5 7.5 10 12.5 15 17.5 20 22 24 30 40 50 60 70 100]';
shape = 'triangle9_s';
for jK = 1:size(s11,1) 
    sp = s11(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_11(1,jK) = imag(Most_unstab)/ut/Rt;
end

% Trapezium
s12 = [5 7.5 10 12.5 15 17.5 20 30 40 50]';
shape = 'trapezium_s';
for jK = 1:size(s12,1) 
    sp = s12(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    points_12(1,jK) = imag(Most_unstab)/ut/Rt;
end

% plot(s1*sqrt(pi/8),points_1,'-o','Color','#E32017','LineWidth',2,'MarkerSize',9);
% plot(s2*sqrt(2/5),points_2,'-+','Color','#003688','LineWidth',2,'MarkerSize',9);
% plot(s3*sqrt(0.5 + sqrt(3)/4),points_3,'-^','Color','#003688','LineWidth',2,'MarkerSize',9);
plot(s4*sqrt(sqrt(3))/2,points_4,'->','LineWidth',2,'Color','#9364CD');
plot(s5/2,points_5,'-v','Color','#E32017','LineWidth',2,'MarkerSize',9);
plot(s6*sqrt(sqrt(3))/2,points_6,'-d','LineWidth',2,'Color','#9364CD');

% plot(s7*sqrt(pi/8),points_7,'--o','Color','#E32017','LineWidth',2,'MarkerSize',9);
% plot(s8*sqrt(2/5),points_8,'--+','Color','#003688','LineWidth',2,'MarkerSize',9);
% plot(s9*sqrt(0.5 + sqrt(3)/4),points_9,'--^','Color','#003688','LineWidth',2,'MarkerSize',9);
plot(s10*sqrt(sqrt(3))/2,points_10,'-->','LineWidth',2,'Color','#9364CD');
plot(s11/2,points_11,'--v','Color','#E32017','LineWidth',2,'MarkerSize',9);
plot(s12*sqrt(sqrt(3))/2,points_12,'--d','LineWidth',2,'Color','#9364CD');

set(gcf,'position',[160 280 800 600])
set(gca,'Xlim',[0 35])
set(gca,'Fontn','Times','FontSize',30,'LineWidth',2)
xlabel('$\ell_g^+$','Interpreter','latex','FontSize',38)
ylabel('$\sigma_{I\,max}^+$','Interpreter','latex','FontSize',38)
% legend('semi-circle','blade','triangle3','triangle6','triangle9','trapezium',...
%     'location','Southeast','FontSize',18)
% title('Maximum amplification vs l_g^+, set at virtual origin')
box on