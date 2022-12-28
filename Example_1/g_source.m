function g = g_source(X,Y,t)
% Original source function $\bar{G}$, to be calculated at time step 't'.
% 'X' and 'Y' are grid data.

x = X; % simplify typing
y = Y;

g = -exp(-t).*(sin(pi*x).*sin(pi*y) + (pi+1).*(x+y) + 1)...
    -(exp(-t)/12).*(2*pi + 2 + pi.*sin(pi.*(x+y)))...
    +(pi^2*exp(-t)./12).*(2+1.5*x+2*y).*sin(pi*x).*sin(pi*y);

g = g(:);

end