%% Programme to solve for Stokes velocity profile in riblets
% Solve A*u = P
% S = riblets grid points
% P = pressure gradient
% A = Laplacian coefficients
% Sc = Sd - S = riblets surface grid points

clear
close all
tic

%% Input data
% Sides of the non-periodic-in-z domain
ly = 0.5; 
lz = 1;
error_hist = [];
n_hist = [];

% Define grid sizes
for n = 5:1:30
    n_hist(end+1) = n;
    ny = n; % number of points y (wall-normal) domain is enclosed
    nz = 2*n; % number of points in z (spanwise) % periodic in z-direction
    dy = ly/(ny-1); % length of sub-intervals in y-axis
    dz = lz/(nz-1); % length of sub-intervals in z-axis
    z = linspace(0,lz,nz);
    y = linspace(0,ly,ny);

    dpdx = -1; % Pressure gradient
    Sx = 0; % Shear at the top, could set to 1 (if normalised, always 1)

    geometry = 2; % 0 = parabola, 1 = triangle, 2 = semi-circle
    savefile = 0;

    %% Build S matrix - Grid Points of riblets
    % S = 1 for points within and on boundary and = 0 elsewhere 
    S = zeros(ny,nz);
    circle = (-sqrt((lz/2)^2-(z-lz/2).^2))+lz/2;
    shape = 'circle';
    for k = 1:nz
        for j = 1:ny
            if y(j)<circle(k)
                S(j,k) = 1;
            else
                S(j,k) = 0;
            end
        end
    end

    %% Build Sd matrix
    % Sd = 1 for points within and on boundary and = 0 elsewhere 
    % riblets now bigger by d = max(dy,dz)
    Sd = zeros(ny,nz);

    circle = (-sqrt((lz/2)^2-(z-lz/2).^2))+lz/2+dy;
    for k = 1:nz
        for j = 1:ny
            if y(j)<circle(k)
                Sd(j,k) = 1;
            else
                Sd(j,k) = 0;
            end
        end
    end

    %% Build Sc matrix - curve of the boundary
    Sc = Sd-S;

    %% Build P matrix - pressure gradient - RHS of Stokes
    P = ones(ny*nz,1);
    P = sparse(P);
    P = P*dpdx;

    %% Build A matrix 
    beta = -2/(dy^2)-2/(dz^2); % coefficient of point p
    alpha = 1/(dy^2); % coefficient of surrounding points
    r1 = beta.*ones(ny,1);
    r2 = alpha.*ones(ny-1,1);
    r3 = alpha.*ones(1,1);
    B = diag(r1,0) + diag(r2,1) + diag(r2,-1) + diag(r3,-ny+1) + diag(r3,ny-1);
    B = sparse(B);
    C = kron(eye(nz),B);
    C = sparse(C);
    gamma = 1/(dz^2);
    r4 = gamma.*ones(ny*nz-ny,1);
    r5 = gamma.*ones(ny,1);
    D = diag(r4,-ny)+diag(r4,+ny)+diag(r5,-ny*nz+ny)+diag(r5,ny*nz-ny);
    D = sparse(D);
    A = C+D;
    A = sparse(A);

    %% Shear at the top condition
    for k = 1:nz
        j = ny; % top of domain
        index = (k-1)*ny+j;
        A(index,:) = 0;
        A(index,index) = 1/dy^2; % should be dy, but P(index) changed as well
        if index~=1
            A(index,index-1) = -1/dy^2;
        end   
        P(index) = Sx/dy;
    end

    %% Solid condition - modify the A matrix
    % Modify A and P matrix for points within and on the boundary of riblets
    % Velocity inside riblets = 0
    for j = 1:ny
        for k = 1:nz
            if S(j,k) == 1
                index = (k-1)*ny+j;
                A(index,:) = 0;
                A(index,index) = 1/dy^2;
                P(index) = 0;
            end
        end
    end

    %% Modify the A matrix - curve of the boundary
    % Modification of the matrix A near the boundary of the solid
    for j=1:ny
        for k=1:nz
            if Sc(j,k) == 1
                y0 = (j-ny)*dy;
                z0 = k*dz ;          
                my = abs(floor(log10(dy)));
                mz = abs(floor(log10(dz)));

            %ys: (j-+2,k);(j-+1,k);(j,k-+2);(j,k-+1) 
                ys=[y0-dy, y0, y0+dy];
                zs=[z0-dz, z0, z0+dz];

                yp0=-(-sqrt((lz/2)^2-(z0-lz/2).^2))+ly;
                zpPlus=sqrt((lz/2)^2-(y0-ly).^2)+lz/2;
                zpMinus=-zpPlus;

                % check closest distance in z direction               
                if abs(z0-zpPlus) < abs(z0 - zpMinus)
                    zp0 = zpPlus;
                else
                    zp0 = zpMinus;
                end

                % Modify ys: (j-+2,k),(j-+1,k),(j,k-+2),(j,k-+1)         
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

    %% No slip at bottom condition
    for k = 1:nz
        j = 1;
        index = (k-1)*ny+j ;
        A(index,:) = 0;
        A(index,index) = 1/dy^2;
        P(index)=0;
    end

    %% Solution
    u = reshape(mldivide(A,P),ny,nz);
    Ub = mean(mean(u));
    % Lws=Ub/s;
    % daspect([1 1 1])
    % axis tight
    % axis fill
    % axis image
    % us = [zeros(1,nz); u];
    us = [u(:,end) u];

    %% Analytical solution for semi-circlular riblets
    % Used to validate results
    za = linspace(0,lz,nz);
    ya = linspace(0,ly,ny);
    nu = 1;
    R = lz/2;
    index = logical(S);
    for j=1:length(ya)
        for k=1:length(za)
            f(j,k)=0.25*(R^2-((ya(j)-ly)^2+(za(k)-lz/2)^2));
    %         fr(j,k)=0.25*(R^2-((ya(j))^2+(za(k))^2));
        end
    end
    Ubmaen=mean(mean(f));
    % Ubmaenr=mean(mean(fr));
    ua=-1/nu*dpdx*f;
    ua(index)=0;
    umax=max(max(ua));
    umaxmax=max(max(u));

    %% Verify
    error = (abs(lz^2/4 - max(u(:))*4))*100/(lz^2/4) % percentage error
    error_hist(end+1) = error;
end
figure
semilogy(n_hist, error_hist)
ylabel('Percentage error')
xlabel('Number or points in y')
toc