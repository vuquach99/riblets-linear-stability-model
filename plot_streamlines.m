% Streamlines plotting script - isocontour of streamfunction for the most amplified mode
clear
close all

nosmod = 512;
Rt = 550;
s = 30;
shape = 'circle1_s';

[D0,D1,D2,D3,D4] = dmat(nosmod);

ymax = [500 100 85 55 40 40 25];
ii = 1;
figure(5)
hold on
for jK = 1:length(s)
%     ii=ii+1;
    sp = s(jK);
    fname = ['Rt' num2str(Rt) '_' shape '_sp' num2str(sp) '_Ny' num2str(nosmod) '.mat']
    load(fname)
    for ia=1:length(lxp)
        maxeig(ia)=Max_unstab(ia)
    end
    [~,i]=max(maxeig);
    lxpi = Most_lxp;
    alp=2*pi*Rt/lxpi;
    maxeigvc=maxeigvc(:,i);

    Nos = size(maxeigvc,1);
    x = 2*lxpi/Rt*[0:200]'/200; % creates x
    yy=Rt*(y+1);
    evx=exp(1i*alp*x');
    x=Rt*x;
    uuu=1i*D1*maxeigvc/alp;
    uuu=real(uuu*evx);
    vvv=D0*maxeigvc;
    vvv=real(vvv*evx);
    phi=1i*D0*maxeigvc/alp;
    phi=real(phi*evx);
%     x=x-49.7;%x=x-x(jj(ii));
    uuu=uuu(yy<ymax(ii),:);
    vvv=vvv(yy<ymax(ii),:);
    phi=phi(yy<ymax(ii),:);
    yy=yy(yy<ymax(ii));

    contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
        [.2:.2:  1],    '-b',   'LineW',2)
    contour(     x(1:2:end),yy(1:2:end),phi(1:2:end,1:2:end) /max(phi(:)), ...
        [-1:.2:-.2],    '-r',   'LineW',2)
     axis equal
    set(gca,'Xlim',[0 lxpi])
    set(gca,'Ylim',[0 15])
end

set(gca,'Fontn','Times','FontSize',34,'linew',2)
xlabel('x^+','FontAngle','italic');
ylabel('y^+','FontAngle','italic');
box on