osmat

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
