% Example 2 data

% domain [0,l_1]x[0,l_2]
l1 = 1;
l2 = 1;

% sensor coordinates (sample)
xs = 0.1*[4  4 16 27 27 39 50 50 62 73 73  85]/10; % cm
ys = 0.1*[13 3  8 13  3  8  13 3  8 13  3   8]/1.5; % cm

% internal measurement domain [x1,x2]x[y1,y2] \subset [0,l_1]x[0,l_2]
x1 = 0.02;
x2 = 0.98;
y1 = 0.2;
y2 = 0.8;

% coefficient function 'C'
CC = @(x,y) 1./(1 + x.^2 + 2*y.^2); 

% exact conductivity
k11 = @(x,y) (3*x.^2 + y.^2 + 1)/2; % k_{11}
k22 = @(x,y) (x.^2 + 3*y.^2 + 1)/2; % k_{22}

% initial temperature
T_sol = @(x,y,t) exp(-pi*t).*cos(pi*t).*(1+x.^2+2*y.^2).*sin(pi*x).*cos(pi*y);

% heat flux functions
F1 = @(y,t) -(pi/2)*exp(-pi*t).*cos(pi*t).*cos(pi*y).*(y.^2+1).*(1+2*y.^2); 
F2 = @(y,t) -pi*exp(-pi*t).*cos(pi*t).*cos(pi*y).*(4+y.^2).*(1+y.^2);
F3 = @(x,t) exp(-pi*t).*cos(pi*t).*sin(pi*x).*(1+x.^2);
F4 = @(x,t) -exp(-pi*t).*cos(pi*t).*sin(pi*x).*(3*x.^2+11);

% source function
G = @(X,Y,t) g_source(X,Y,t);

% heat transfer functions
H1 = @(y) ones(size(y));
H2 = @(y) ones(size(y));
H3 = @(x) ones(size(x));
H4 = @(x) ones(size(x));




















