% Adequates data to our working set [0,l_1]x[0,l_2].
% run 'data.m' first.

c = @(x,y) CC(l1*x,l2*y);

k = @(x,y) Kond(l1*x,l2*y);

T0 = 22*ones(size(X));
T0 = T0(:);

f1 = @(y,t) F1(l2*y,t); 
f2 = @(y,t) F2(l2*y,t);
f3 = @(x,t) F3(l1*x,t); 
f4 = @(x,t) F4(l1*x,t);

g = @(X,Y,t) G(l1*X,l2*Y,t);

h1 = @(y) H1(l2*y);
h2 = @(y) H2(l2*y);
h3 = @(x) H3(l1*x);
h4 = @(x) H4(l1*x);