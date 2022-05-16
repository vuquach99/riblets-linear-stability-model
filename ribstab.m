% Extracts relevant results and processes data

%% Initialises matrices to store eigenvectors and eigenvalues
eigvecs = zeros(nosmod+1,nosmod+1,nx); % (eigenvector's 2 dimensions, wavenumber)
eigvals = zeros(nosmod+1,nx); % (eigenvalue, wavenumber)

%% Stores most amplified eigenvalue of each wavenumber
Max_unstab = zeros(1,nx);
for ia = 1:nx
    alp = alp0(ia);
    d = turchan(alp,D0,D1,D2,D3,D4,U,nut,Kvp,Kup,Kvs,Kus);
    [eigvecs(:,:,ia),eigvals(:,ia)] = iord2(d);
    
    % Finds omega_i >= 10 and makes those eigenvalues NaN
    Max_imag1 = find(imag(eigvals(:,ia)) >= 10);
    for i = 1:length(Max_imag1)
        eigvals(Max_imag1(i),ia) = NaN;
    end
    
    % Stores eigenvalues with -6 <= omega_i & 0 < omega_r < 100 in Unstab
    Max_imag1 = find(imag(eigvals(:,ia))>=-6);
    Max_real1 = zeros(1,length(Max_imag1));
    for in = 1:length(Max_imag1)
        index = Max_imag1(in);
        Max_real1(in) = real(eigvals(index,ia));
        if Max_real1(in) < 100 && Max_real1(in) > 0
            Unstab(in,ia) = eigvals(index,ia); % Unstab = matrix of wavenumber's eigenvalues (columns)
        else
            Unstab(in,ia) = -1000-1000i; % Other eigenvalues are pushed to this value
        end
    end
    % If Unstab(:,ia) is empty, populate it with zeros
    if isempty(Max_imag1)
        Unstab = [Unstab(:,:) zeros(size(Unstab,1),1)];
    end
    
    % Then stores the most amplified eigenvalue of this wavenumber in Max_unstab
    [~,Max_ind] = max(imag(Unstab(:,ia)));
    Max_unstab(ia) = Unstab(Max_ind,ia);
end

% Stores the most amplified eigenvalue of all wavenumbers in Most_unstab
[~,Most_ind] = max(imag(Max_unstab));
Most_unstab = Max_unstab(Most_ind);
Most_lxp = lxp(Most_ind);

%% Post-processing
maxeigvl=[];
maxeigvc=[];
Vmaxeig =[];
for i=1:nx
    lxpi=lxp(i);
    alp=2*pi*Rt/lxpi;
    ii=(1:size(eigvals,1));
    eevl=squeeze(eigvals(ii,i));
    eevc=squeeze(eigvecs(:,ii,i));
    if isempty(eevl)
        maxeigvl(i)=0;
        maxeigvc(:,i)=0*(1:size(y,1));
        Vmaxeig (i)=0;
    elseif isempty(find(max(imag(eevl))))
        indmax=1;
        maxeigvl(i)=eevl(indmax);
        maxeigvc(:,i)=eevc(:,indmax);
        Vmaxeig(i)=real(eevl(indmax))/alp/ut;
    else
        indmax=find(max(imag(eevl)));
        maxeigvl(i)=eevl(indmax);
        maxeigvc(:,i)=eevc(:,indmax);
        Vmaxeig(i)=real(eevl(indmax))/alp/ut;
    end
end
