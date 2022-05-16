function d = turchan(alp,D0,D1,D2,D3,D4,U,nut,Kvp,Kup,Kvs,Kus)

% Creates Orr-Sommerfeld matrices.
% Uses Chebyshev pseudospectral discretisation for channel flow.

% INPUT:
% alp = alpha (streamwise wave number)
% Dn = nth derivative matrix

% OUTPUT:
% A,B = Orr-Sommerfeld matrices
% d/dt(B*v)=-i*A*v ==> d=inv(B)*A;

ak2 = alp^2;
Nos = size(U,1);

Nutpp = nut(:,3);
Nutp = nut(:,2);
Nut = nut(:,1);

u = U(:,1);
uy = U(:,2);
uyy = U(:,3);
unos = ones(1,length(u));

% Sets up Orr-Sommerfeld matrices
B =  D2 - ak2*D0;

A = alp*(u*unos).*B;
A = A - alp*(uyy*unos).*D0;
A = A + 1i*(Nut*unos).*(D4 - 2*ak2*D2 + (ak2^2)*D0);
A = A + 2i*(Nutp*unos).*(D3 - ak2*D1);
A = A + 1i*(Nutpp*unos).*(D2 + ak2*D0);

% Boundary condition 1
Bt = D1(1,:); % top layer
Bb = D1(end,:); % bottom layer

At = alp*u(1)*D1(1,:) + (-alp*uy(1) - 1i*Nut(1)/Kvp)*D0(1,:);
At = At + 1i*Nut(1)*Kvs/Kvp*D2(1,:);
At = At + 1i*Nut(1)*(D3(1,:) - ak2*D1(1,:));
At = At + 1i*Nutp(1)*(D2(1,:) + ak2*D0(1,:));

Ab = alp*u(end)*D1(end,:) + (-alp*uy(end) + 1i*Nut(end)/Kvp)*D0(end,:);
Ab = Ab - 1i*Nut(end)*Kvs/Kvp*D2(end,:); % +?
Ab = Ab + 1i*Nut(end)*(D3(end,:) - ak2*D1(end,:));
Ab = Ab + 1i*Nutp(end)*(D2(end,:) + ak2*D0(end,:));

% Boundary condition 2
Atu = Kup/Kvp*(-Kvs*D2(1,:) - D0(1,:)) - D1(1,:) + Kus*D2(1,:); % -Kup/Kvp... -Kus?
Abu = -Kup/Kvp*(Kvs*D2(end,:) + D0(end,:)) - D1(end,:) + Kus*D2(end,:); % +Kup/Kvp... -D0?
Btu = Atu;
Bbu = Abu;

% Pushes the spurious eigenvalues somewhere far away
er = -20000i;
B = [Bt; Btu; B(3:Nos-2,:); Bbu; Bb];
A = [At; er*Atu; A(3:Nos-2,:); er*Abu; Ab];

% Returns d
d = B\A;
end
