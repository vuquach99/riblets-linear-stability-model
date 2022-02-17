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
jj=[94];
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
for ii=1:1
  set(gca,'Fontn','Times','FontSize',24,'linew',2)
  set(gca,'XTick',[0:xtcksp(ii):800])
  set(gca,'YTick',[0:ytcksp(ii):200])
  xlabel('x^+','Ver','bo','Hor','r','FontAngle','italic');
  ylabel(' y^+  ','Ver','c','Hor','c','FontAngle','italic');
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
hold on
  for jK=[7]
    ii=ii+1
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
    x=x-49.7;%x=x-x(jj(ii));
    uuu=uuu(yy<ymax(ii),:);
    vvv=vvv(yy<ymax(ii),:);
    phi=phi(yy<ymax(ii),:);
    yy=yy(yy<ymax(ii));
    contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
        [.2:.2:  1],    '-b',   'LineW',2)
    contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
        [-1:.2:-.2],    '-r',   'LineW',2)
    set(gca,'Xlim',[0 lxpi])
    set(gca,'Ylim',[0 .3582*lxpi])
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print -depsc2 new_pert_streamlines1_COLOR
