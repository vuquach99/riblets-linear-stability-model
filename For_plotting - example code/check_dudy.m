nosmod=1536;
Kwp=10;
Rtt=550;
NK=size(Kwp,1);
NR=size(Rtt,1);
zi=sqrt(-1);

eval(['load basefuns2_' num2str(nosmod) ';'])
nm=round(N/2);

jj=5;%[1 30 94];
ymax=40;
xtcksp=40;
ytcksp=15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1),clf,hold on
set(gcf,'position',[1 440 1250 560])
set(gcf,'PaperPositionMode','auto')

ii=1
  set(gca,'Fontn','Times','FontSize',20)
  set(gca,'XTick',[0:xtcksp(ii):800])
  set(gca,'YTick',[0:ytcksp(ii):200])
  xlabel('x^+','Ver','bo','Hor','r','FontAngle','italic');
  ylabel(' y^+  ','Ver','c','Hor','c','FontAngle','italic');
  axis equal
  box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jRe=1
  Rt=Rtt(jRe)
  jK=1
  ii=1
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

   dudy=zeros(size(uuu,1)-1,size(uuu,2))
   dy=(yy(1:end-1)-yy(2:end))
   for i=1:size(dudy,2)
     dudy(:,i)=(uuu(1:end-1,i)-uuu(2:end,i))./(yy(1:end-1)-yy(2:end));
   end

    mylincontour(x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end)'/max(phi(:)), ...
        [.2:.2:  1],'','--k',10,'linew',2)
    contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
        [-1:.2:-.2],    '-k',   'LineW',2)
    set(gca,'Xlim',[0 lxpi])
    set(gca,'Ylim',[0 .3582*lxpi])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%print -depsc2 aaaaaaaaaaaa
