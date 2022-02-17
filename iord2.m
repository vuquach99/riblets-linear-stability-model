function [xs,es]=iord2(d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   compute ordered eigenvalues of d              %
%   d==input matrix                               %
%   es==ordered eigenvalues                       %
%   xs==ordered eigenvectors                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [v,e]=eig(d);
  e=diag(e);
  [eimag,is]=sort(-imag(e));
  xs=v(:,is);
  es=e(is);
