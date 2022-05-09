function [xs,es] = iord2(d)
% Computes ordered eigenvalues of d

% INPUT:
% d = input matrix

% OUTPUT:
% es = ordered eigenvalues
% xs = ordered eigenvectors

[v,e] = eig(d); 
% v = matrix with eigenvectors as columns
% e = diagonal matrix of eigenvalues
e = diag(e);
[~,is] = sort(imag(e),'descend');
xs = v(:,is);
es = e(is);
