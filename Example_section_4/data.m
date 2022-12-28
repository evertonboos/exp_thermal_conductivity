% Example data

% domain [0,l_1]x[0,l_2]
l1 = 0.1; % m
l2 = 0.015; % m

% sensor coordinates
xs = 0.001*[4,27,50,73,4,27,50,73,16,39,62,85]; % x-coordinates (sensor)
ys = 0.001*[3,3,3,3,13,13,13,13,8,8,8,8]; % y-coordinates (sensor)

%xs = 0.001*[4  4 16 27 27 39 50 50 62 73 73  85]; % cm x-coordinates (sensor)
%ys = 0.001*[13 3  8 13  3  8  13 3  8 13  3   8]; % cm y-coordinates (sensor)

% internal measurement domain [x1,x2]x[y1,y2] \subset [0,l_1]x[0,l_2]
x1 = min(xs);
x2 = max(xs);
y1 = min(ys);
y2 = max(ys);

% final time
tf = 60; % seconds

% coefficient function 'C'
CC = @(x,y) (1.18e5*44.5)*ones(size(x)); % m2/s

% exact conductivity
Kond = @(x,y) 44.5*ones(size(x)); % k; % W/m°C

% heat flux functions
f0 = 21;
F1 = @(y,t) f0*ones(size(y)); 
F2 = @(y,t) f0*ones(size(y)); 
F3 = @(x,t) f0*ones(size(x)); 
F4 = @(x,t) f0*ones(size(x)); 

% source function
G = @(X,Y,t) g_source(X,Y,t);

% heat transfer functions
%ho = 0.002866; % cal/s cm2 Celsius
ho = 120; % W/m2°C
H1 = @(y) ho*ones(size(y));
H2 = @(y) ho*ones(size(y));
H3 = @(x) ho*ones(size(x));
H4 = @(x) ho*ones(size(x));





