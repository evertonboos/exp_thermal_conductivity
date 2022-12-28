% Example data

% domain [0,l_1]x[0,l_2]
l1 = 1;
l2 = 1;

% internal measurement domain [x1,x2]x[y1,y2] \subset [0,l_1]x[0,l_2]
x1 = 0.02;
x2 = 0.98;
y1 = 0.2;
y2 = 0.8;

% coefficient function 'C'
CC = @(x,y) ones(size(x)); 

% exact conductivity
k11 = @(x,y) (1 + x + y)/12; % k_{11}
k22 = @(x,y) (1 + 0.5*x + y)/12; % k_{22}

% k22 = @(x,y) (1 + x + y)/12; % k_{11}
% k11 = @(x,y) (1 + 0.5*x + y)/12; % k_{22}

% exact solution
T_sol = @(x,y,t) exp(-t).*(sin(pi*x).*sin(pi*y) + (pi+1).*(x+y) + 1);

% heat flux functions
F1 = @(y,t) T_sol(0,y,t) - ((1+y)/12).*exp(-t).*(pi.*sin(pi*y)+pi+1); 
F2 = @(y,t) T_sol(1,y,t) + ((2+y)/12).*exp(-t).*(-pi.*sin(pi*y)+pi+1);
F3 = @(x,t) T_sol(x,0,t) - ((1+0.5*x)/12).*exp(-t).*(pi.*sin(pi*x)+pi+1);
F4 = @(x,t) T_sol(x,1,t) + ((2+0.5*x)/12).*exp(-t).*(-pi.*sin(pi*x)+pi+1);

% F3 = @(y,t) T_sol(0,y,t) - ((1+y)/12).*exp(-t).*(pi.*sin(pi*y)+pi+1); 
% F4 = @(y,t) T_sol(1,y,t) + ((2+y)/12).*exp(-t).*(-pi.*sin(pi*y)+pi+1);
% F1 = @(x,t) T_sol(x,0,t) - ((1+0.5*x)/12).*exp(-t).*(pi.*sin(pi*x)+pi+1);
% F2 = @(x,t) T_sol(x,1,t) + ((2+0.5*x)/12).*exp(-t).*(-pi.*sin(pi*x)+pi+1);

% source function
G = @(X,Y,t) g_source(X,Y,t);

% heat transfer functions
H1 = @(y) ones(size(y));
H2 = @(y) ones(size(y));
H3 = @(x) ones(size(x));
H4 = @(x) ones(size(x));

% H3 = @(y) ones(size(y));
% H4 = @(y) ones(size(y));
% H1 = @(x) ones(size(x));
% H2 = @(x) ones(size(x));





