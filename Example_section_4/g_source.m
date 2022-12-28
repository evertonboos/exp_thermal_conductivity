function g = g_source(X,Y,t)
% Original source function $\bar{G}$, to be calculated at time step 't'.
% 'X' and 'Y' are grid data.
% Normal distribution source

%% Internal parameters
Wx = 1.91; % J
Wy = 0.91;
Wz = 3.25;
l1 = 0.1; % m
l2 = 0.015; % m
%Ar = 1.5e-6; % m2
Ar = 1.5e-7;
%abar = 0.0001; % m
vf = 1.7e-3; % m/s in [0,l1]
tf = 60; % s (source application until the end of the piece)
q0 = (Wx+Wy+Wz)/(Ar*tf); % W/m^2, or J/m^2 s

%% Source evaluation at mesh [X,Y] and time = t
x = X(:,1); % 'x' Chebyshev's points
y = Y(1,:); % 'y' Chebyshev's points
xx = X(:); % ordenate x grid
yy = Y(:); % ordenate y grid
n = length(x)-1; % number of Chebyshev's points
g = 0.0*ones(n+1,n+1);
G = g;
y0 = 0.015; % to center the source on the side of the piece

% source evaluation
v = vf*t;
x0 = vf*t;

if x0 == 0 || x0 > l1
    g = zeros((n+1)^2,1);
else
    g = q0*(36/pi)*exp(-(3*(yy-y0/2)./(y0/2)).^2).*exp(-(3*(xx-x0)./x0).^2);
end

g = g(:);

end




