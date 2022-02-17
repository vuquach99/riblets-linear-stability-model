	function [D0,D1,D2,D3,D4,y]=Dmat(N);

%
% Function to create differentiation matrices
%
% y	= wall-normal coordinate
% N	= number of modes
% D0	= zeroth derivative matrix
% D1	= first  derivative matrix
% D2	= second derivative matrix
%
	y = cos(pi*(0:N)/(N))';

% initialize

	num = round(abs(N));

% create D0

	D0 = [];
	vec = (0:1:num)';
	for j=0:1:num
	   D0 = [D0 cos(j*pi*vec/num)];
	end 

% create higher order derivative matrices

	lv = length(vec);
	D1 = [zeros(lv,1) D0(:,1)     4*D0(:,2)  ];
	D2 = [zeros(lv,1) zeros(lv,1) 4*D0(:,1)  ];
    D3 = [zeros(lv,1) zeros(lv,1) zeros(lv,1)]; 
	D4 = [zeros(lv,1) zeros(lv,1) zeros(lv,1)];
	for j=3:num;
	   D1 = [D1 2*j*D0(:,j)+j*D1(:,j-1)/(j-2)];
	   D2 = [D2 2*j*D1(:,j)+j*D2(:,j-1)/(j-2)];
       D3 = [D3 2*j*D2(:,j)+j*D3(:,j-1)/(j-2)];
	   D4 = [D4 2*j*D3(:,j)+j*D4(:,j-1)/(j-2)];
	end
