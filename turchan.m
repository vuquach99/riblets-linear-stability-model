function d = turchan(alp,D0,D1,D2,D3,D4,U,Rt,iKw,nut,iSw,Re,lhq,lsq,L)
% will change to: (alp,D0,D1,D2,D3,D4,U,nut,iKw,iSw,lhq,lsq)

% Creates Orr-Sommerfeld matrices.
% Uses Chebyshev pseudospectral discretisation for channel flow.

% INPUT:
% alp = alpha (streamwise wave number)
% Dn = nth derivative matrix

% OUTPUT:
% A,B = Orr-Sommerfeld matrices
% d/dt(B*v)=-zi*A*v    ==>   d=inv(B)*A;

% a=iKw/iSw.*D0(1,:)
% some parameters:
ak2 = alp^2;
Nos = size(U,1);

nutpp=nut(:,3); nutp=nut(:,2); nut=nut(:,1);
% uses mean velocity and turbulent viscosity from turprof
% uses analytical expression from Reynolds and Tiedermann

% set up Orr-Sommerfeld matrix

Uy=U(:,2); Uyy=U(:,3); U=U(:,1);
unos=ones(1,length(U));

% d/dt(B*v)=-zi*(A*v) ==> d=inv(B)*A;  %
B11 =  D2 - ak2*D0;

A11 =  alp*(U*unos).*B11;                                  %%% conveccion
A11 =  A11 - alp*(Uyy*unos).*D0;                           %%% U''v
A11 =  A11 + 1i*(nut*unos).*(D4 - 2*ak2*D2 + (ak2^2)*D0);
A11 =  A11 + 1i*(nutp*unos).*( 2*D3 - 2*ak2*D1);  %%% viscous, nutp
A11 =  A11 + 1i*(nutpp*unos).*( D2 + ak2*D0);

% boundary conditions
% v
% A11t= - iKw * D0(1,:) - iSw * D1(1,:);
% A11b= + iKw * D0(end,:) - iSw * D1(end,:);
% 
% B11t= A11t ; %%%  d/dt*(dv/dy)=(i*alp*Uy-nu/Kw)*v
% B11b= A11b;

B11t= D1(1  ,:) ; %%%  d/dt*(dv/dy)=(i*alp*Uy-nu/Kw)*v
B11b= D1(end,:) ; %%%  d/dt*(dv/dy)=(i*alp*Uy+nu/Kw)*v
% 
A11t= (-alp*Uy(1) - 1i*iKw)*D0(1,:) ;
A11t= A11t + alp*U(1)*D1(1,:); % U =/= 0 added
% % A11t= A11t -i*1/Re*(D3(1,:)-ak2*D1(1,:));
A11t= A11t + 1i*nut(1)*(D3(1,:)-ak2*D1(1,:));
A11t= A11t + 1i*nutp(1)*(D2(1,:)+ak2*D0(1,:));
%With shear 
A11t= A11t + 1i*D2(1,:)*L;

A11b= (-alp*Uy(end) + 1i*iKw)*D0(end,:) ;
A11b= A11b + alp*U(end)*D1(end,:);
% % A11b= A11b -i*1/Re*(D3(end,:)-ak2*D1(end,:));
A11b= A11b + 1i*nut(end)*(D3(end,:)-ak2*D1(end,:));
A11b= A11b + 1i*nutp(end)*(D2(end,:)+ak2*D0(end,:));
% With shear
A11b= A11b + 1i*D2(end,:)*L;

% A11t= (-alp*Uy(2)-i*iKw)*D0(2,:) ;
% A11t= A11t + alp*U(2)*D1(2,:); 
% A11t= A11t +i*nut(2)*(D3(2,:)-ak2*D1(2,:));
% 
% A11b= (-alp*Uy(end-1)+i*iKw)*D0(end-1,:) ;
% A11b= A11b + alp*U(end-1)*D1(end-1,:);
% A11b= A11b +i*nut(end-1)*(D3(end-1,:)-ak2*D1(end-1,:));

% u
% PRESSURE
% A11tu= - iKw * D0(1,:) - iSw * D1(1,:);
% A11bu= + iKw * D0(end,:) - iSw * D1(end,:);
% 
% B11tu= A11tu ; %%%  d/dt*(dv/dy)=(i*alp*Uy-nu/Kw)*v
% B11bu= A11bu; %%%  d/dt*(dv/dy)=(i*alp*Uy+nu/Kw)*v

% SHEAR
A11tu = -(-D2(1,:)*lhq-D0(1,:))*iKw/iSw;
A11tu = A11tu - D1(1,:) + D2(1,:)*lsq;

A11bu = (D2(end,:)*lhq-D0(end,:))*iKw/iSw;
A11bu = A11bu - D1(end,:) - D2(end,:)*lsq;

B11tu= A11tu ; %%%  d/dt*(dv/dy)=(i*alp*Uy-nu/Kw)*v
B11bu= A11bu ;
% % % 
% % % B11tu= D1(1,:) ; %%%  d/dt*(dv/dy)=(i*alp*Uy-nu/Kw)*v
% % % B11bu= D1(end,:) ; %%%  d/dt*(dv/dy)=(i*alp*Uy+nu/Kw)*v
% % 
% A11tu= -alp*Uy(1)*D0(1,:);
% A11tu= A11tu + i*iSw*D1(1,:) ;
% A11tu= A11tu + alp*U(1)*D1(1,:); 
% A11tu= A11tu + i*nut(1)*(D3(1,:)-ak2*D1(1,:));
% 
% A11bu= -alp*Uy(end)*D0(end,:);
% A11bu= A11bu + i*iSw*D1(end,:) ;
% A11bu= A11bu + alp*U(end)*D1(end,:);
% A11bu= A11bu + i*nut(end)*(D3(end,:)-ak2*D1(end,:));
% % 
% B11 = [ B11t ; B11(2:Nos-1,:) ; B11b ];
% A11 = [ A11t ; A11(2:Nos-1,:) ; A11b ];

er  = -20000i; 
B11 = [ B11t ; B11tu ; B11(3:Nos-2,:) ; B11bu ; B11b ];
A11 = [ A11t ; er*A11tu; A11(3:Nos-2,:); er*A11bu; A11b ];

% % Only shear driven flow
% B11 = [D0(1,:)-D2(1,:) ; D1(1,:)-D2(1,:); B11(3:Nos-2,:) ; D1(Nos,:)-D2(Nos,:); D0(Nos,:)-D2(Nos,:) ];
% A11 = [er*(D0(1,:)-D2(1,:)) ; er*(D1(1,:)-D2(1,:)); A11(3:Nos-2,:); er*(D1(Nos,:)-D2(Nos,:)); er*(D0(Nos,:)-D2(Nos,:)) ];

% 
% B11 = [ B11tu ; B11t ; B11(3:Nos-2,:) ; B11b ; B11bu ];
% A11 = [ A11tu ; A11t; A11(3:Nos-2,:); A11b; A11bu ];

% % du/dy = 0;
% B11 = [ B11t ; D2(1,:); B11(3:Nos-2,:) ; D2(Nos,:); B11b ];
% A11 = [ A11t ; er*D2(1,:); A11(3:Nos-2,:); er*D2(Nos,:); A11b];

d=B11\A11;
