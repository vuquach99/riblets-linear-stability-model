clear
clc
q = 0.35;
s=[20]';%1 5 10 14 16 [1 2 5 10 11 14 16 18 20 25 30 40 50]'; %[5 10 20 30 40 50]';
Lwp=q*((pi/256)^(1/3))*s;
Lsp=q*((1/24)^(1/2))*s; 
%% Triangular riblets
% F1=0.004448936237952;
% F2=0.020842780639986;
%% Trapezoidal riblets
% F1=0.007178406796948;
% F2=0.027606404326489;
%% Blade
% F1=0.0083;
% F2=0.0229;
%% Permeabilities expression 
% Lwp=F1^(1/3)*s;
% Lsp=F2^(1/2)*s;


G1=0.1269;
G2=0.2507;
Lhq=q*s*(G1^(1/2));
Lsq=q*s*G2;

Rtt=[550]';
% Sweep through a range of wavelengths
nosmod=256;%256%
nx=50;
lxpmin=5;
lxpmax=10000;
% For investigating a single wavelength value
% nx= 1; % number of wavelengths
% lxpmin=56;
% lxpmax=57;
    
kapa=0.426;
Aint =25.4;
eddyfrac=1;

[D0,D1,D2,D3,D4]=Dmat_inviscid(nosmod); % Chebychev polynomials

lxp=log(lxpmin):(log(lxpmax)-log(lxpmin))/(nx-1):log(lxpmax);
lxp=fliplr(exp(lxp));
zi=sqrt(-1);

NK=size(Lwp,1);
NR=size(Rtt,1);
for jR=1:NR
  Rt=Rtt(jR)
  [y,nut,uc,Re]=turprof_generic(nosmod,Aint,kapa,eddyfrac,Rt); % Cess turbulent velocity profile inside channel
  utau=Rt/Re;
  alp0=2*pi*Rt./lxp;
  for jK=1:NK

    Lwpp=Lwp(jK); % wall-normal coeff. in wall units
    Kw=(Lwpp/Rt)^3; % wall normal coeff. in outer units cubed (Lw3)
    iKw=1/Kw/Re;
    
    Lspp=Lsp(jK); % streamwise coeff. in wall units
    Sw=(Lspp/Rt)^2; % streamwise coeff. in outer units squared (Ls2)
    iSw=1/Sw/Re;
    
    % Coefficients for shear driven flow
    lhq=(Lhq(jK)/Rt)^2; 
    lsq=Lsq(jK)/Rt;
    L=lhq/Kw/Re;

    fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat'];
    ribstab_inviscid_noD_noProf
    save(fname,'Rt','Max_unstab','y','Lwpp','nosmod','utau','lxp','maxeigvc','maxeigvl','eigvals')
  end
end

n = length(Lwp);
figure
for r=1:length(Lwp)
    Lwpp=Lwp(r);
    fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat']
    load(fname)
    a{r}= imag(eigvals)/utau/Rt>0.2;
    b{r} = imag(eigvals(a{r}));
%     maxim=max(imag(eigvals));
%     maxim=max(maxim)/utau/Rt;
%     plot(real(eigvals)/utau/Rt,imag(eigvals)/utau/Rt,'.', 'Color', [(r-1)*1/n', 0, 1-(r-1)*1/n],'MarkerSize', 10)
     plot(real(eigvals)/utau/Rt,imag(eigvals)/utau/Rt,'.', 'Color', [(r-1)*1/n', 0, 1-(r-1)*1/n],'MarkerSize', 10)
    hold on
end

xlabel('\omega^+_R')
ylabel('\omega^+_I')
title('Orr-Sommerfeld spectrum')
box on
set(gca,'Fontn','Times','FontSize',14,'LineWidth',2)

%% Most amplified mode spectrum comparison between press and shear+press driven flow
% figure
hold on
mrk{1}='.'; mrk{2}='x'; mrk{3}='k-s'; mrk{4}='k-^'; %'k-o' 'k-*'
clr(1,:)='k'; clr(2,:)='k'; clr(3,:)='k'; clr(4,:)='k';
line{1} = 'k-'; line{2} = 'k--'; line{3} = 'k:';
    marker=mrk{1};
    for r=1:length(Lwp)
        Lwpp=Lwp(r);
        fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat'];
        load(fname)
        plot(real(Max_unstab)/utau/Rt,imag(Max_unstab)/utau/Rt,marker, 'Color', [(r-1)*1/n', 0, 1-(r-1)*1/n],'MarkerSize', 22)
        hold on
%         Max_unstab_all(r,c)=Max_unstab;
    end
% Difference=(abs(imag(Max_unstab_all(:,1))-imag(Max_unstab_all(:,2))))
% indexes=find(max(Difference));
% Difference(indexes)=NaN
% figure
% Diff=[0.2409 0.5886 0.9532 2.4845]
% Rtt=[180 550 1000 2000];
% plot(Rtt, Diff)
% legend=ftypes;
xlabel('\omega^+_R')
ylabel('\omega^+_I')
% title('Most amplified mode comparison')
box on
set(gca,'Fontn','Times','FontSize',14,'LineWidth',2)