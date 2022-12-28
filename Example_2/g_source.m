function g = g_source(X,Y,t)
% Original source function $\bar{G}$, to be calculated at time step 't'.
% 'X' and 'Y' are grid data.

x = X; % simplify typing
y = Y;

% repeating terms
K11 = (3*x.^2 + y.^2 + 1)/2;
K22 = (x.^2 + 3*y.^2 + 1)/2;
Ct = 1 + x.^2 + 2*y.^2;
expt = exp(-pi*t).*cos(pi*t);
sx = sin(pi*x);
sy = sin(pi*y);
cx = cos(pi*x);
cy = cos(pi*y);
%exp_x = exp(-pi*t).*cos(pi*t).*sin(pi*x);
%exp_y = exp(-pi*t).*cos(pi*t).*cos(pi*y);

%g = -pi.*( cos(pi*t) + sin(pi*t) ).*exp(-pi*t).*sin(pi*x).*cos(pi*y) ...
%    -K11.*( 4*pi*x.*cos(pi*x)-(pi^2*Ct-2).*sin(pi*x) ).*exp_y ...
%    -K22.*( -(pi^2*Ct-4).*cos(pi*y) -8*pi*y.*sin(pi*y) ).*exp_x...
%    -3*x.*( pi*Ct.*cos(pi*x)+2*x.*sin(pi*x) ).*exp_y...
%    -3*y.*( 4*y.*cos(pi*y)-pi*Ct.*sin(pi*y) ).*exp_x;

g = -pi*exp(-pi*t).*( cos(pi*t) + sin(pi*t) ).*sx.*cy...
    -3*x.*expt.*cy.*( pi*cx.*Ct + 2*x.*sx )...
    -K11.*expt.*cy.*( -pi^2*sx.*Ct + 4*pi*x.*cx + 2*sx )...
    -3*y.*expt.*sx.*( -pi*sy.*Ct + 4*y.*cy )...
    -K22.*expt.*sx.*( -pi^2*cy.*Ct - 8*pi*y.*sy + 4*cy );

%g = flip(g);
g = g(:);

end