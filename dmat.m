function [D0,D1,D2,D3,D4,y]=dmat(N)
% Creates differentiation matrices

% INPUT: 
% N	= number of modes

% OUTPUT:
% y	= wall-normal coordinate
% Dn = nth derivative matrix

y = cos(pi*(0:N)/(N))';

% initialises
num = round(abs(N));

% creates D0
D0 = [];
vec = (0:1:num)';
for j = 0:1:num
   D0 = [D0 cos(j*pi*vec/num)];
end 

% creates higher order derivative matrices
lv = length(vec);
D1 = [zeros(lv,1) D0(:,1)     4*D0(:,2)  ];
D2 = [zeros(lv,1) zeros(lv,1) 4*D0(:,1)  ];
D3 = [zeros(lv,1) zeros(lv,1) zeros(lv,1)]; 
D4 = [zeros(lv,1) zeros(lv,1) zeros(lv,1)];
for j = 3:num
   D1 = [D1 2*j*D0(:,j)+j*D1(:,j-1)/(j-2)];
   D2 = [D2 2*j*D1(:,j)+j*D2(:,j-1)/(j-2)];
   D3 = [D3 2*j*D2(:,j)+j*D3(:,j-1)/(j-2)];
   D4 = [D4 2*j*D3(:,j)+j*D4(:,j-1)/(j-2)];
end
