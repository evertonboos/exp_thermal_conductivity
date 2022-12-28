%% CHEBY1  compute D = differentiation matrix, x = Chebyshev grid
  function [D,y] = cheby1(N,a,b);
  %a =1; b = 1;
   if N==0, D=0; x=1; return, end
   x = cos(pi*(0:N)/N)'; y=0.5*((b-a)*x+(b+a));
  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
   X = repmat(y,1,N+1);
  dX = X-X';
   D1  = (c*(1./c)')./(dX+(eye(N+1)));   % off-diagonal entries
   D1  = D1 - diag(sum(D1'));          % diagonal entries
   y=flipud(y);                       % changing points from left to right
   D = flipud(fliplr(D1));
