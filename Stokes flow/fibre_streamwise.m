clear
close all
tic
%
%% Input data

% Define fibre in unit cell
nfibres = 3;
nzi = 36;
nyi = 36;
a = 0.05;%rad(iii);%0.10; % ellipsoidal fibre, major semiaxis
b = a; % ellipsoidal fibre, minor semiaxis
theta = 0; % ellipsoidal fibre, inclination angle
% Define domain
ly = nfibres+1;
lz = 1;
% Define grid sizes
ny_p = nyi*ly; % number of points (y) 
ny = ny_p+1; % number of points (y)
nz = nzi; % number of points (z)
dy = ly/(ny_p); % length of sub-intervals in y-axis
dz = lz/nz; % length of sub-intervals in z-axis
yaxis = (-ly:dy:0);
zaxis = (0:dz:lz);

% Pressure gradient
dpdx = 0;
% Shear at the top
Sx = 1;

%% Generate fibres data set - ellipsoidal fibers:

yc = -ly+0.5*ones(1,nfibres); % y-centre of the fibre
zc = lz/2; % z-centre of the fibre
for i = 2:nfibres
    yc(i) = yc(i-1)+1;
end
rd=max(dy,dz);

%% Build S matrix - [Shape of the body dependent]
% SS = 1 for points within and on boundary of fibre and = 0 elsewhere 
% SS stored in S_ for each fibre

SS = zeros(ny,nz);
for i = 1:nfibres
    for j=1:ny
        for k=1:nz
            y=(j-(ny_p+1))*dy;
            z=k*dz;    
            SS(j,k)=(1-sign(round( (((z-zc)*sin(theta)-(y-yc(i))*cos(theta))^2)/a^2+(((z-zc)*cos(theta)+(y-yc(i))*sin(theta))^2)/b^2-1    ,10)))/2+isinf(1/(sign(round( (((z-zc)*sin(theta)-(y-yc(i))*cos(theta))^2)/a^2+(((z-zc)*cos(theta)+(y-yc(i))*sin(theta))^2)/b^2-1  ,10))))/2;
        end
    end 
    S_{i} = SS;
    S_{i} = sparse(S_{i});
end
clear SS

%% Build Sd matrix (solid with r'=r+d) - [Shape of the body dependent]
% SSd = 1 for points within and on boundary of fibre with radius r' = r+d and = 0 elsewhere 
% where r = radius of fibre and d = max(dy,dz)
% Sd_ stored in S_ for each fibre

SSd = zeros(ny,nz);
for i = 1:nfibres
    for j=1:ny
        for k=1:nz
            y=(j-(ny_p+1))*dy;
            z=k*dz;
            SSd(j,k)=(1-sign(round( (((z-zc)*sin(theta)-(y-yc(i))*cos(theta))^2)/(a+rd)^2+(((z-zc)*cos(theta)+(y-yc(i))*sin(theta))^2)/(b+rd)^2-1    ,10)))/2+isinf(1/(sign(round( (((z-zc)*sin(theta)-(y-yc(i))*cos(theta))^2)/(a+rd)^2+(((z-zc)*cos(theta)+(y-yc(i))*sin(theta))^2)/(b+rd)^2-1  ,10))))/2;
        end
    end 
    Sd_{i} = SSd;
    Sd_{i} = sparse(Sd_{i});
    Sc_{i} = Sd_{i} - S_{i};
end
clear Sd_ SSd


%% Build P matrix - pressure gradient 

P = ones(ny*nz,1);
P = sparse(P);
P = dpdx*P;

%% Build A matrix 

beta = (-2/(dy^2)-2/(dz^2)); 
alpha = 1/(dy^2); 
r1 = beta.*ones(ny,1);
r2 = alpha.*ones(ny-1,1);
r3 = alpha.*ones(1,1);
B = diag(r1,0) + diag(r2,1) + diag(r2,-1)+ diag(r3,-ny+1) + diag(r3,ny-1) ;
B = sparse(B);
C = kron(eye(nz),B);
C = sparse(C);
gamma = 1/(dz^2);
r4 = gamma.*ones(ny*nz-ny,1);
r5 = gamma.*ones(ny,1);
D=diag(r4,-ny)+diag(r4,+ny)+diag(r5,-ny*nz+ny)+diag(r5,ny*nz-ny);
D = sparse(D);
A=C+D;
A=sparse(A);

%% Solid condition - modify the A matrix
% Modify A and P matrix for points within and on the boundary of each fibre

for i = 1:nfibres
    S = cell2mat(S_(i));
    for j = 1:ny
        for k = 1:nz
            if S(j,k)==1
                index =(k-1)*ny+j ;
                A(index,:) = 0;
                A(index,index) = 1;
                P(index) =0;
            end
        end
    end
end
clear S

%% Modify the A2 matrix
% Modification of the matrix A near the boundary of the solid

for i = 1:nfibres
    S = cell2mat(S_(i));
    Sc = cell2mat(Sc_(i));
    ycc = yc(i);
    for j=1:ny
        for k=1:nz
            if Sc(j,k)==1
                y0=(j-(ny_p+1))*dy;
                z0=k*dz ;          
                my=abs(floor( log10(dy)));
                mz=abs(floor( log10(dz)));
            
            %ys: (j-+2,k);(j-+1,k);(j,k-+2);(j,k-+1) 
                ys=[y0-dy, y0, y0+dy];
                zs=[z0-dz, z0, z0+dz];         
                           
            %Find the points of intersection: (j-+1,k),(j,k-+1) 
                pz1=(((sin(theta))^2)/a^2+((cos(theta))^2)/b^2);
                pz2=(-2*zc*(((sin(theta))^2)/a^2+((cos(theta))^2)/b^2)-2*(y0-ycc)*sin(theta)*cos(theta)*(1/a^2-1/b^2));
                pz3=(zc^2*(((sin(theta))^2)/a^2+((cos(theta))^2)/b^2)+(y0-ycc)*2*zc*sin(theta)*cos(theta)*(1/a^2-1/b^2)+(y0-ycc)^2*(((cos(theta))^2)/a^2+((sin(theta))^2)/b^2)-1);
                pz=[pz1 pz2 pz3];
                Rz=roots(pz);
                Rz=isreal(Rz)*Rz+(1-isreal(Rz))*(z0-2*dz);
                if abs(z0-Rz(1))<abs(z0-Rz(2))
                    zp0=Rz(1);
                else
                    zp0=Rz(2);
                end

                py1=(((sin(theta))^2)/b^2+((cos(theta))^2)/a^2);          
                py2=(-2*ycc*(((sin(theta))^2)/b^2+((cos(theta))^2)/a^2)-2*(z0-zc)*sin(theta)*cos(theta)*(1/a^2-1/b^2));         
                py3=(ycc^2*(((sin(theta))^2)/b^2+((cos(theta))^2)/a^2)+(z0-zc)*2*ycc*sin(theta)*cos(theta)*(1/a^2-1/b^2)+(z0-zc)^2*(((cos(theta))^2)/b^2+((sin(theta))^2)/a^2)-1);
                py=[py1 py2 py3];
                Ry=roots(py);
                Ry=isreal(Ry)*Ry+(1-isreal(Ry))*(y0-2*dy);               
                if abs(y0-Ry(1))<abs(y0-Ry(2))
                    yp0=Ry(1);
                else
                    yp0=Ry(2);
                end      
            
                %Modify ys: (j-+2,k),(j-+1,k),(j,k-+2),(j,k-+1)         
                deltay=round(abs(y0-yp0),10);
                deltays=y0-yp0;
                deltaz=round(abs(z0-zp0),10);
                deltazs=z0-zp0;
            
                if deltaz<dz
                    if round(deltaz,mz+2)~=0
                        zs(2-sign(deltazs))=z0-sign(deltazs)*deltaz;
                    end
                end           
                if deltay<dy
                    if round(deltay,my+2)~=0 
                        ys(2-sign(deltays))=y0-sign(deltays)*deltay;                     
                    end
                end 
            
                bj=zeros(1,3);
                bk=zeros(1,3);          
                bj(1)=2/((ys(2)-ys(1))*(ys(3)-ys(1)));
                bj(2)=-2/((ys(3)-ys(2))*(ys(2)-ys(1))); 
                bj(3)=2/((ys(3)-ys(2))*(ys(3)-ys(1)));
                bk(1)=2/((zs(2)-zs(1))*(zs(3)-zs(1)));
                bk(2)=-2/((zs(3)-zs(2))*(zs(2)-zs(1))); 
                bk(3)=2/((zs(3)-zs(2))*(zs(3)-zs(1)));  
            
                index =(k-1)*ny+j;
                A(index,index)= bj(2) + bk(2);        
                if abs(deltay)<dy           
                    if j==1
                    A(index,index+ny-1)= bj(1);
                    A(index,index+1)= 0;
                  
                    elseif j==ny
                    A(index,index-ny+1)= bj(3);
                    A(index,index-1)= 0;

                    else
                    A(index,index-sign(deltays))=0;       
                    A(index,index+sign(deltays))= bj(2+sign(deltays));

                    end
                end                                        
                if abs(deltaz)<dz 
                    if k==1
                    A(index,index+(nz-1)*ny)= bk(1);
                    A(index,index+ny)= bk(3);
                    elseif k==nz
                    A(index,index-ny)= bk(3);
                    A(index,index-(nz-1)*ny)= bk(1);
                    else
                    A(index,index-sign(deltazs)*ny)= bk(2-sign(deltazs));
                    A(index,index+sign(deltazs)*ny)= bk(2+sign(deltazs));
                    end
                end
            end
        end
    end
end
%% Boundary condition
% Shear at the top
for k=1:nz
    j=ny;
    index =(k-1)*ny+j ;
    A(index,:) = 0;
    A(index,index) = 1/dy;
    A(index,index-1) = -1/dy;
    P(index)=Sx;
end
% No slip at the bottom
for k = 1:nz  %controllare qui
    j=1;
    index =(k-1)*ny+j ;
    A(index,:) =0;
    A(index,index) =1;
%     A(index,index) = -2/dz^2+1/dy^2;
%     A(index,index+1) = -2/dy^2;
%     A(index,index+2)= 1/dy^2;
    P(index)=0;
end

%% Solution
u = reshape(mldivide(A,P),ny,nz);

% us = [zeros(1,nz); u];
us = [u(:,end) u];

%% Verify
for k=1:nz
    Sxtot(k)=(u(ny,k)-u(ny-1,k))/dy;
end
Sxtotm=mean(Sxtot);

%% average in z
ums=mean(u,2);
% ums = [0; um];
ums = full(ums);

% %% Figure

figure(1)
contourf(zaxis,yaxis, full(us))
hold on
t=0:0.01:2*pi; 
for i = 1:nfibres
    zpatch{i} = b*cos(t)*cos(theta)-a*sin(t)*sin(theta)+zc;
    ypatch{i} = b*cos(t)*sin(theta)+a*sin(t)*cos(theta)+yc(i);
    patch(zpatch{i},ypatch{i},'w')
end
axis image
grid on
set(gca,'layer','top')
title('u')
colorbar
ylabel('y')
xlabel('z')
%
screen = nfibres;
screenbot = ny_p*(screen-1)/ly+1;
screentop = ny_p*screen/ly+1;
figure(11)
hold on
pcolor(zaxis,yaxis(screenbot:screentop),full(us(screenbot:screentop,:)))
shading interp
patch(zpatch{screen},ypatch{screen},'w')
axis image
colorbar
ylabel('y')
xlabel('z')
title('u')
set(gca,'layer','top')
%
figure(2)
hold on
plot( ums, yaxis,'k','Linewidth',1.25)
yline(yc(end)+a,'--')
axis image
% grid on
box on
ylabel('$y$','Interpreter','Latex')
xlabel('$\bar{u}$','Interpreter','Latex')
set(gca,'FontName','times','FontSize',15)
% title('')
jinterface = find(yaxis>=(yc(end)+a),1);
p = polyfit(yaxis(jinterface:end),ums(jinterface:end),2);
plot(polyval(p,yaxis(jinterface-5:end)),yaxis(jinterface-5:end),'r')
% tau_interface = 2*p(1)*yaxis(jinterface)+p(2)
tau_interface = (ums(jinterface+1)-ums(jinterface))/(yaxis(jinterface+1)-yaxis(jinterface))



%%
toc

