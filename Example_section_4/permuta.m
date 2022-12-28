 
%
% calcula permutacao que transforma enumeracao de 
% de baixo para acima e de esquerad a direita, em
% de esquerda a dieita e de baixo para cima 
%  per: satisfaz:  per*per = I; per'*per =I=per*per';
% vale per*kron(I,A)*per = kron(A,I);
function [x,per] = permuta(n)
 n = n+1; 
 x = 1:n*n;
 a = reshape(x,n,n); a = a';
 p = a(:); Id= eye(n);
 I  = eye(n*n);
 per = sparse(I(p,:));
 
 %A = randn(4);
 %Dx = kron(Id,A);
 %Dy = kron(A,Id);
 %[D,u] = ncheb(3,1,0);