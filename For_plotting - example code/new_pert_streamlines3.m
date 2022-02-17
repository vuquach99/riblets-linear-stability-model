nosmod=1536;
%Lwp=[1.5 3 6 12 24 48 96]';
Kwp=[2 4 6 8 10 12 14]';
Rtt=[185 550 950 2000]';
NK=size(Kwp,1);
NR=size(Rtt,1);
zi=sqrt(-1);

eval(['load basefuns2_' num2str(nosmod) ';'])
nm=round(N/2);

%jj=[1 52 30 5 37 94];
jj=[1 30 94];
ymax=[100 85 55 40 40 25];
xtcksp=[40 40 40 40 40 40];
ytcksp=[15 15 15 15 15 15];

pos(1,:)=[0.1300    0.5838    0.2134    0.3412];
pos(2,:)=[0.4108    0.5838    0.2134    0.3412];
pos(3,:)=[0.6916    0.5838    0.2134    0.3412];
pos(4,:)=[0.1300    0.3000    0.2134    0.3412];
pos(5,:)=[0.4108    0.3000    0.2134    0.3412];
pos(6,:)=[0.6916    0.3000    0.2134    0.3412];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1),clf
set(gcf,'position',[1 440 1250 560])
set(gcf,'PaperPositionMode','auto')

for ii=1:3
  subplot(2,3,ii)
  set(gca,'Fontn','Times','FontSize',22,'linew',2.2)
  set(gca,'XTick',[0:xtcksp(ii):800])
  set(gca,'YTick',[0:ytcksp(ii):200])
%  xlabel('x^+','Ver','bo','Hor','r','FontAngle','italic');
%  ylabel(' y^+  ','Ver','c','Hor','c','FontAngle','italic');
end
subplot(2,3,1),hold on
xlabel('$x^+$','interpreter','latex','Pos',[21 -1 0]);
ylabel('$y^+$','interpreter','latex','Pos',[-1.5 8 0]);
subplot(2,3,2),hold on
xlabel('$x^+$','interpreter','latex','Pos',[25 -1.4 0]);
ylabel('$y^+$','interpreter','latex','Pos',[-2 8 0]);
subplot(2,3,3),hold on
xlabel('$x^+$','interpreter','latex','Pos',[25 -1.4 0]);
ylabel('$y^+$','interpreter','latex','Pos',[-2 8 0]);
for ii=1:3
  subplot(2,3,ii)
  axis equal
  box on
end
%subplot(2,3,1)
%xlabel('            x^+','Ver','bo','Hor','l','FontAngle','italic');
%subplot(2,3,2)
%xlabel('            x^+','Ver','bo','Hor','l','FontAngle','italic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jRe=2
  Rt=Rtt(jRe)
  ii=0
  for jK=[2 4 7]
    ii=ii+1
    subplot(2,3,ii)
    Kwpp=Kwp(jK);
    fname=['ribstab_Rt' num2str(Rt) '_Kw' num2str(Kwpp) '_Ny' num2str(nosmod) '.mat'];
    load(fname)
    [Amax,i]=max(imag(maxeigvl));
    lxpi=lxp(i)
    alp=2*pi*Rt/lxpi;
    maxeigvc=maxeigvc(:,i);
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
    x=x-x(jj(ii));
    uuu=uuu(yy<ymax(ii),:);
    vvv=vvv(yy<ymax(ii),:);
    phi=phi(yy<ymax(ii),:);
    yy=yy(yy<ymax(ii));
    mylincontour(x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end)'/max(phi(:)), ...
        [.2:.2:  1],'','--k',10,'linew',2)
    contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
        [-1:.2:-.2],    '-k',   'LineW',2)
    set(gca,'Xlim',[0 lxpi])
    set(gca,'Ylim',[0 .3582*lxpi])
    set(gca,'Position',pos(ii,:))
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print -depsc2 new_pert_streamlines3
