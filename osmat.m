% Computes Orr-Sommerfeld / Squire matrix for 3D Poiseuille Flow and to compute energy matrix.

eigvals = zeros(         nosmod+1,nx);
eigvecs = zeros(nosmod+1,nosmod+1,nx);

% Sweeps through wavenumbers
for ia = 1:nx
    alp = alp0(ia);
    d = turchan(alp,D0,D1,D2,D3,D4,U,nut,Kvp,Kup,Kvs,Kus);
    [eigvecs(:,:,ia),eigvals(:,ia)]=iord2(d);

Max_imag1=find(imag(eigvals(:,ia))>=10);

for i=1:length(Max_imag1)
 eigvals(Max_imag1(i),ia)=NaN;
end

Max_imag1=find(imag(eigvals(:,ia))>=-6);
Max_real1=zeros(1,length(Max_imag1));
for in=1:length(Max_imag1)
    index=Max_imag1(in);
    Max_real1(in)=real(eigvals(index,ia));
    if Max_real1(in)<100 && Max_real1(in)>0
        Unstab(in,ia)=eigvals(index,ia);
    else
        Unstab(in,ia)=-1000-1000i;
    end
end
[value, Max_ind]=(max(imag(Unstab)));
Max_unstab=Unstab(Max_ind);

end