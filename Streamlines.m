% script to plot isocontour of streamfunction for the most amplified mode

nosmod=512; %

[D0,D1,D2,D3,D4]=Dmat_inviscid(nosmod);

% s=[5];
s=[20]';%[1 2 6 10 11 14 16 50 60 70 80]';%[1 2 6 8 10 11 12 14 16 18 20 25 30 40 50]';
flowtype='viscous';
F1=0.012273817101988 % (pi/256);
F2=(1/24);
%% Triangular riblets
% F1=0.004448936237952;
% F2=0.020842780639986;
%% Trapezoidal riblets
% F1=0.007178406796948;
% F2=0.027606404326489;
%% Blade
% F1=0.0083;
% F2=0.0229;
%% Permeabilities
Lwp=F1^(1/3)*s;
Lsp=F2^(1/2)*s;

Rtt=550; %

zi=sqrt(-1);
ftypes={'BCtest'}

% eval(['load basefuns2_' num2str(nosmod) ';'])
% nm=round(N/2);

%jj=[1 52 30 5 37 94];
ymax=[100 85 55 40 40 25];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure,clf
for ii=1:1
  set(gca,'Fontn','Times','FontSize',24,'linew',2)
%   set(gca,'XTick',[0:xtcksp(ii):800])
%   set(gca,'YTick',[0:ytcksp(ii):200])
  xlabel('x^+','Ver','bo','Hor','r','FontAngle','italic');
  ylabel(' y^+  ','Ver','c','Hor','c','FontAngle','italic');
  axis equal
  box on
end

jRe=1
  Rt=Rtt(jRe)
%   ii=0
hold on
flowtype=ftypes{1}
  for jK=1:length(Lwp)
%     ii=ii+1
    Lwpp=Lwp(jK);
    fname=['ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat']
% fname=['results\ribstab_Rt' num2str(Rt) '_Lw' num2str(Lwpp) '_Ny' num2str(nosmod) '.mat']
    load(fname)
    for k=1:length(lxp)
        maxeig(k)=max(imag(eigvals(:,k)));
    end
    [Amax,i]=max(maxeig);
%     [Amax,i]=max(imag(eigvals));% 
% [Amax,i]=max(imag(maxeigvl));
% [Amax,i]=max(imag(mineigvl));
    lxpi=lxp(i)
    alp=2*pi*Rt/lxpi;
    maxeigvc=maxeigvc(:,i);
% mineigvc=mineigvc(:,i);
    Nos = size(maxeigvc,1);
    x=2*lxpi/Rt*[0:200]'/200;
    yy=Rt*(y+1);
    evx=exp(zi*alp*x');
    x=Rt*x;
    uuu=zi*D1*maxeigvc/alp;
    uuu=real(uuu*evx);
    vvv=D0*maxeigvc;
    vvv=real(vvv*evx);
    phi=zi*D0*maxeigvc/alp;
    phi=real(phi*evx);
    x=x-49.7;%x=x-x(jj(ii));
    uuu=uuu(yy<ymax(ii),:);
    vvv=vvv(yy<ymax(ii),:);
    phi=phi(yy<ymax(ii),:);
    yy=yy(yy<ymax(ii));
% %%     Poiseulle Test
% d2pdx2=nut(:,1).*uuu*D2*alp;
% RHS=d2pdx2+nut(:,1)*(D2*uuu+alp^2*uuu);
% LHS=(-zi*maxeigvl+zi*alp*uc(:,1))+zi*alp*uuu*uc(:,2);
%%
    contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
        [.2:.2:  1],    '-b',   'LineW',2)
    contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
        [-1:.2:-.2],    '-r',   'LineW',2)
    set(gca,'Xlim',[0 lxpi])
    set(gca,'Ylim',[0 .3582*lxpi])
  end
  
% title([num2str(ftypes{1}) ' driven flow streamlines'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print -depsc2 new_pert_streamlines1_COLOR

% 
% nosmod=512;
% %Lwp=[1.5 3 6 12 24 48 96]';
% %Lwp=[2 4 6 8 10 12 14]';
% s=[2 5 10 15 20]'; % riblet spacing in wall units [2 4 6 8 10 12 14]';
% Lwp=(pi/256)*s; %.^3;
% Rtt=[550] % 550 950 2000]';
% NK=size(Lwp,1);
% NR=size(Rtt,1);
% zi=sqrt(-1);
% 
% % eval(['load basefuns2_' num2str(nosmod) ';'])
% % nm=round(N/2);
% nm=round(nosmod/2);
% 
% %jj=[1 52 30 5 37 94];
% jj=[94];
% ymax=[100 85 55 40 40 25];
% xtcksp=[40 40 40 40 40 40];
% ytcksp=[15 15 15 15 15 15];
% 
% pos(1,:)=[0.1300    0.5838    0.2134    0.3412];
% pos(2,:)=[0.4108    0.5838    0.2134    0.3412];
% pos(3,:)=[0.6916    0.5838    0.2134    0.3412];
% pos(4,:)=[0.1300    0.3000    0.2134    0.3412];
% pos(5,:)=[0.4108    0.3000    0.2134    0.3412];
% pos(6,:)=[0.6916    0.3000    0.2134    0.3412];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(1),clf
% for ii=1:1
%     set(gca,'Fontn','Times','FontSize',24,'linew',2)
%     set(gca,'XTick',[0:xtcksp(ii):800])
%     set(gca,'YTick',[0:ytcksp(ii):200])
%     xlabel('x^+','Ver','bo','Hor','r','FontAngle','italic');
%     ylabel(' y^+  ','Ver','c','Hor','c','FontAngle','italic');
%     axis equal
%     box on
% end
% %subplot(2,3,1)
% %xlabel('            x^+','Ver','bo','Hor','l','FontAngle','italic');
% %subplot(2,3,2)
% %xlabel('            x^+','Ver','bo','Hor','l','FontAngle','italic');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% jRe=length(Rtt)
% Rt=Rtt(jRe)
% ii=0
% hold on
% for jK=1:length(Lwp)
%     ii=ii+1
%     Lwpp=Lwp(jK);
%     for i=1:nx
%         fname=['results/nosmod_' num2str(nosmod) '_Rt_' num2str(Rt) '_Lwp_' num2str(Lwpp) '.mat'];
%         load(fname)
%         [Amax,i]=max(imag(eigvals(i,:))); %.*Re./Rt.^2;
%         lxpi=lxp(i)
%         alp=2*pi*Rt/lxpi;
%         maxeigvcs=max(maxeigvc(i,:));
%     end
%     Nos = size(maxeigvcs,1);
%     x=2*lxpi/Rt*[0:512]'/512;
%     yy=Rt*(y+1);
%     evx=exp(zi*alp*x);
%     x=Rt*x;
%     uuu=zi*D1*maxeigvcs/alp;
%     uuu=real(uuu*evx);
%     vvv=D0*maxeigvcs;
%     vvv=real(vvv*evx);
%     phi=zi*D0*maxeigvcs/alp;
%     phi=real(phi*evx);
%     x=x-49.7;%x=x-x(jj(ii));
%     uuu=uuu(yy<ymax(ii),:);
%     vvv=vvv(yy<ymax(ii),:);
%     phi=phi(yy<ymax(ii),:);
%     yy=yy(yy<ymax(ii));
%     contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
%         [.2:.2:  1],    '-b',   'LineW',2)
%     contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
%         [-1:.2:-.2],    '-r',   'LineW',2)
%     set(gca,'Xlim',[0 lxpi])
%     set(gca,'Ylim',[0 .3582*lxpi])
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% print -depsc2 new_pert_streamlines1_COLOR
% 
% 
% %         maxeigvc=[];
% %         ii=(1:size(eigvals,1));
% %         eevl=squeeze(eigvals(ii,i));
% %         eevc=squeeze(eigvecs(:,ii,i));
% %         if isempty(eevl)
% %             maxeigvl(i)=0;
% %             maxeigvc(:,i)=0*(1:size(y,1));
% %         elseif isempty(find(max(imag(eevl))))
% %             indmax=1;
% %             maxeigvl(i)=eevl(indmax);
% %             maxeigvc(:,i)=eevc(:,indmax);
% %         else
% %             indmax=find(max(imag(eevl)));
% %             maxeigvl(i)=eevl(indmax);
% %             maxeigvc(:,i)=eevc(:,indmax);
% %         end
