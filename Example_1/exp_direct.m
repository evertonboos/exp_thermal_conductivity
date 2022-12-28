% Direct problem of finding T given k(x,y) and modelled data.

clear; close all;

%% Parameters and grid details
n = 15; % Chebyshev's points 
N = 13; % number of equally spaced interval on time discretization
tf = 1; % final time

t0 = linspace(0,tf,N+1);
tt0 = t0;
et = 0;
%t0 = [0, t];
alpha = 1e-1;
beta = 1-alpha;

% plot positioning
aa = (1920-1366)/1920*0.5; % notebook -> monitor display
bb = (1080-768)/1080*0.5;
%posi = [0 0 1 1]; % full screen (any monitor)
posi = [aa bb 1366/1920 768/1080]; % centered, notebook size

[D,x] = cheby1(n,0,1); % differentiation matrix 
y = x; % grid
[~,P] = permuta(n); % permutation matrix
hat_D = [alpha*D(:,1), D(:,2:end-1), alpha*D(:,end)]; % \hat{D}
d0 = D(:,1); 
dn = D(:,n+1);
d0_e1 = sparse(n+1,n+1); 
d0_e1(:,1) = d0; % product d_0*e_{1}^T
dn_en = sparse(n+1,n+1); 
dn_en(:,end) = dn; % product d_n*e_{n+1}^T

[X,Y] = meshgrid(x,y);
X = X';
Y = Y';

% data allocation
mydata = struct('cheby', {n,N,t0,tt0,et,alpha,beta,D,x,y,P,hat_D,d0,dn,d0_e1,dn_en,X,Y});
run('data.m') % original data
run('data_adequation.m') % data adequation to [0,1]x[0,1]

K11 = k11(X,Y);
K22 = k22(X,Y);
K = [K11(:); K22(:)];

%% Generate temperature
TK = ort_generate_T(K,T0,mydata);

figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;
Err = zeros(1,N+1);
for i = 1:N+1
    subplot(1,2,1)
    Te = T(X,Y,t0(i));
    Err(i) = norm(TK(:,i)-Te(:))/norm(Te(:));
    %Te = g(X,Y,t0(i)); Te = reshape(Te,n+1,n+1);
    surf(X,Y,Te)
    title('Exact')
    axis square
    zlim([0, 10])
    %zlim([-.5, 0.5])
    subplot(1,2,2)
    surf(X,Y,reshape(TK(:,i),n+1,n+1))
    error(i) = norm(Te(:) - TK(:,i));
    %surf(X,Y,reshape(g(X,Y,t0(i)),n+1,n+1))
    %zlim([-5,7]), view(140,33)
    title(['Approximation, time = ' num2str(t0(i)) ' '])
    xlabel('x')
    axis square
    zlim([0, 10])
    %zlim([21.5 23])
    %subplot(1,2,2), contour(X2,Y2,T2)
    %axis square
    pause(0.3)
end























