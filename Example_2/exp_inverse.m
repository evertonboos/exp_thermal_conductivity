% Inverse problem of finding k(x,y) given temperature values and modelled
% data.

clear; close all;

sol = zeros(11,1); % solution data save, sol = [n,N,tf,delta,alpha,err1,err_int1,err2,err_int2,iter,errT,DP]'

%% PARAMETERS AND GRID DETAILS
n = 15; % Chebyshev's points 
N = 10; % number of equally spaced intervals on time discretization
tf = 1; % final time
delta = 0.015; % noise level
data_source = 0; % 'exact' to use exact temperature; 'numerical' to use numerical temperature
points = 'full'; % 'full' to use full Cheby's mesh; 'inner' to use only inner domain points; 'selection' to a especial selection of points
chosenNorm = 2; % 2, 'inf'
over = 1; % over repetition on the measurements, 1 for 'full', 4 for 'inner' and 13 for 'selection'
tau = 1.1; % DP parameter
kmax = 50; % maximum number of iterations
alpha = 0; % influence of conductivity on the boundary

sol(1:5) = [n,N,tf,delta,alpha]; % save data

% plot positioning
aa = (1920-1366)/1920*0.5; % notebook -> monitor display
bb = (1080-768)/1080*0.5;
%posi = [0 0 1 1]; % full screen (any monitor)
posi = [aa bb 1366/1920 768/1080]; % centered, notebook size

% time vectors
t0 = linspace(0,tf,N+1); % given time steps
et = 0; % quantity of points to introduce between every given time step
tt0 = linspace(0,tf,(N+1)+N*et); % add extra points inside every interval in t0

beta = 1-alpha;
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

% Chebyshev's grid
[X,Y] = ndgrid(x,y);
xx = X(:);
yy = Y(:);

% data allocation
mydata = struct('cheby', {n,N,t0,tt0,et,alpha,beta,D,x,y,P,hat_D,d0,dn,d0_e1,dn_en,X,Y});
run('data.m') % original data
run('data_adequation.m')% data adequation to [0,1]x[0,1]

%% POINTS SELECTION
% internal routine to identify points to be used during minimization
% process, i.e, select correct indexes of the mesh to be used. Until now,
% possible options are:
% 1) 'full': select full Chebyshev's mesh, to standard minimization procedure
% 2) 'inner' (MS1): select only points inside a given inner rectangular domain (defined in 'data.m')
% 3) 'selection' (MS2): special selection of points (only as example).

if strcmp(points,'full') % FULL DOMAIN SELECTION
    ell_inner = 1:(n+1)^2;
elseif strcmp(points,'inner') % INNER DOMAIN ROUTINE
    xi = zeros(n+1,1);
    yi = zeros(n+1,1);
    for i = 1:n+1
        if l1*x(i) >= x1 && l1*x(i) <= x2
            xi(i) = 1;
        end
        if l2*y(i) >= y1 && l2*y(i) <= y2
            yi(i) = 1;
        end
    end

    ell_inner = zeros((n+1)^2,1); % to select inner Jacobian lines
    for i = 1:n+1
        for j = 1:n+1
            if xi(j) == 1 && yi(i) == 1 
                ell = (j-1) + (i-1)*(n+1) + 1;
                ell_inner(ell) = ell;
            end
        end
    end
    ell_inner = nonzeros(ell_inner); 
elseif strcmp(points,'selection') % SPECIAL SELECTION
    % take n = 15.
    ell_inner = [51,54,57,60,84,87,90,93,117,120,123,126,147,150,153,156,180,183,186,189]';
end
ell_inner = sort(ell_inner); % to reorder the points
ell = length(ell_inner); % number of points selected

% plot domain points
figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;
rectangle('position',[0 0 l1 l2])
hold on
rectangle('position',[x1 y1 x2-x1 y2-y1],'facecolor',[.8 .8 .8],'edgecolor',[.8 .8 .8])
plot(l1*X,l2*Y,'k*')
%plot(xs,ys,'bp','markersize',7)
plot(l1*xx(ell_inner),l2*yy(ell_inner),'r*')
%axis equal
xlim([0-0.1 l1+0.1])
ylim([0-0.1 l2+0.1])
axis square
title('Measurement locations')
hold off

% exact conductivity
Ke11 = k11(X,Y);
Ke22 = k22(X,Y);
Ke = [Ke11(:); Ke22(:)]; 
Kee = Ke;

%% TEMPERATURE MEASUREMENTS AND NOISE
% compute temperatures
Uex = zeros((n+1)^2,N);
for i = 1:N
    TT = T(X,Y,t0(i+1)); % calculate exact temperature at t_{i+1}
    Uex(:,i) = TT(:); % store exact temperature
end
U = Uex;
disp(['DETAILS: Approach for "' points '" case, using exact temperature with ' num2str(100*delta) '% of noise'])
disp('Iterations:')

% noise introduction
rng('default')
Uell = U(ell_inner,:); % exact data in restricted points
Uell = Uell(:); % vectorized given data
Uell = Uell*ones(1,over); % repeat given data 'over' times
%nU = norm(U(:),chosenNorm);
nU = norm(Uell(:),chosenNorm);
%noise = randn((n+1)^2*N,1);
noise = randn(ell*N*over,1); % noise
noise = noise/norm(noise,chosenNorm); % normalized noise vector
errU = delta*nU; 
%Up = U + delta*nU*reshape(noise,(n+1)^2,N); % perturbed U, ||Up-U|| = delta*||U|| (as vector)
Up = Uell(:) + delta*nU*noise; % perturbed U, ||Up-U|| = delta*||U|| (as vector)
%Up = Up(ell_inner,:);
%Up = Up(:);

% stopping criteria setting
DP = tau*errU; % DP criterion

%% SOLUTION FOR 'K' CALCULATED BY LEVENBERG-MARQUARDT'S METHOD
figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;

gamma = 0.7; % LMM constant, \gamma \in (0,1)
sigma = 1e-2; % Armijo's constant, \sigma \in (0,1)

% first iterate K
%Ks = (1)*ones(2*(n+1)^2,1);
v = randn(2*(n+1)^2,1); v = v/max(abs(v));
Ks = 1 + 0.015*v;

% LMM scaling matrix
L = get_l(n+1,2); % first or second order discrete derivative
LL = [kron(eye(n+1),L); kron(L,eye(n+1))];
%D = speye(2*(n+1)^2); % identity matrix, standard case
D = kron(eye(2),LL); % scaling with first or second order derivative 

% storage errors and residuals
res = zeros(kmax,1); 
rel = res;
ERR = zeros(kmax,4); % [err(K11), int_err(K11), err(K22), int_err(K22)]

% ITERATIONS
Kc = zeros(2*(n+1)^2,kmax+1); % storage computed solutions
Kc(:,1) = Ks;
flag = 0;
i = 1; % counter
while (i <= kmax && flag == 0)
    JJ = ort_Jacobian(Ks,mydata); % Jacobian at time steps t0
    J1 = sparse(ell*N,2*(n+1)^2); % inner Jacobian
    for j = 1:N-1
        int = (j-1)*ell+1:j*ell;
        J1(int,:) = JJ(ell_inner,:,j); 
    end
    
    J = zeros(over*ell*N,2*(n+1)^2); % Jacobian for overdetermination
    for j = 1:over
        int = (j-1)*N*ell+1:j*N*ell;
        J(int,:) = J1;
    end
    
    TK_ant = ort_generate_T(Ks,T0,mydata); % solution for K
    %TK_ant = TK_ant(:,2:end);
    TK_ant = TK_ant(ell_inner,1 + (1:N)*(et+1)); % select time steps and measure points
    UK_ant = TK_ant(:); % vectorized solution for K
    UK_ant = UK_ant*ones(1,over); UK_ant = UK_ant(:);
    F_ant = norm(UK_ant-Up);
    mu = F_ant^2;
    
    % compute direction d_k (LMM)
    d = ( J'*J + mu*(D'*D) )\ ( -J'*(UK_ant-Up) ); % direction d_k
    
    TK_pos = ort_generate_T(Ks+d,T0,mydata);
    TK_end = TK_pos(:,end);
    %TK_pos = TK_pos(:,2:end);
    TK_pos = TK_pos(ell_inner,1 + (1:N)*(et+1)); % select time steps and measure points
    UK_pos = TK_pos(:);
    UK_pos = UK_pos*ones(1,over); UK_pos = UK_pos(:);
    F_pos = norm(UK_pos-Up);

    %% New Iterate Definition (LMM with line search)
    if F_pos <= gamma*F_ant % new iterate (without line search)
        Ks = Ks + d;
        step = 1; 
    elseif flag == 0 % LINE SEARCH
        m = 1;
        TK = ort_generate_T(Ks + gamma^m*d,T0,mydata);
        %TK = TK(:,2:end);
        TK = TK(ell_inner,1 + (1:N)*(et+1)); % select time steps and measure points
        UK = TK(:);
        UK = UK*ones(1,over); UK = UK(:);
        F = norm(UK-Up);
        gradphi_d = (J'*(UK_ant-Up))'*d; 
        while 0.5*F^2 - 0.5*F_ant^2 > sigma*gamma^m*gradphi_d % Armijo
            if gamma^m < 1e-10
                disp('Armijo failed, aborting...')
                flag = 1; % abortion flag
                break
            end
            m = m+1;
            TK = ort_generate_T(Ks + gamma^m*d,T0,mydata);
            TK_end = TK(:,end);
            %TK = TK(:,2:end);
            TK = TK(ell_inner,1 + (1:N)*(et+1)); % select time steps and measure points
            UK = TK(:);
            UK = UK*ones(1,over); UK = UK(:);
            F = norm(UK-Up);
        end
        if flag == 0 % correct update if Armijo worked
            Ks = Ks + gamma^m*d; % new iterate
            F_ant = F;
            UK_pos = UK;
            step = gamma^m; % Armijo's step size
        end
    end
    
    Kc(:,i+1) = Ks; % storage, if needed
    normE(i) = norm(Ks-Kee);
    normK(i) = norm(Ks);
    
    %% Stopping criteria
    % relative residual test
    if i >= 2
        %rel(i) = abs(res(i) - res(i-1))/res(i-1);
        rel(i) = norm(Kc(:,i)-Kc(:,i-1))/norm(Kc(:,i));
        if rel(i) <= 5e-4
            disp('-------------------------------------------------------------')
            disp(['Relative residual criteria reached! Iterations: ' num2str(i) ' '])
            flag = 3;
        end
    end
    
    % DP test
    res(i) = norm(UK_pos-Up,chosenNorm); % residual
    if res(i) <= DP
        disp('-------------------------------------------------------------')
        disp(['DP reached! Iterations: ' num2str(i) ' '])
        flag = 2;
    end
    
    %% Displays
    for jj = 1:2
        if jj == 1 % K11
            K = Ks(1:(n+1)^2);
            K = reshape(K,n+1,n+1); K = K(:);
            Ke = Ke11;
            ii = 0; % counter for plotting
        elseif jj == 2 % K22
            K = Ks((n+1)^2+1:end);
            K = reshape(K,n+1,n+1); K = K(:);
            Ke = Ke22;
            ii = 5; % counter for plotting
        end
        err = norm(Ke(:) - K)/norm(Ke(:));
        Ke_int = Ke(2:end-1,2:end-1);
        K_full = reshape(K,n+1,n+1);
        K_int = K_full(2:end-1,2:end-1);
        err_int = norm(Ke_int(:)-K_int(:))/norm(K_int(:));
        ERR(i,2*(jj-1)+1) = err;
        ERR(i,2*(jj-1)+2) = err_int;
        subplot(2,5,ii+2)
        contour(X,Y,reshape(Ke,n+1,n+1),20);
        axis square
        title('exact (contour)')
        subplot(2,5,ii+3)
        contour(X,Y,reshape(K,n+1,n+1),20);
        axis square
        title(['approx., res = ' num2str(res(i)) ', DP = ' num2str(DP) ' '])
        subplot(2,5,ii+4)
        surf(X,Y,reshape(K,n+1,n+1)), zlim([0, 3])
        axis square
        title(['approx., err = ' num2str(err) ', int. err = ' num2str(err_int) ''])
        xlabel('x')
        subplot(2,5,ii+1)
        surf(X,Y,reshape(Ke,n+1,n+1)), zlim([0, 3])
        axis square
        title(['exact, iteration = ' num2str(i) ''])
    end
    T_end = reshape(Uex(:,end),n+1,n+1); % exact temperature at final time
    TK_end = reshape(TK_end,n+1,n+1); % computed temperature at final time
    subplot(2,5,5) % plot computed temperature at final time
    surf(X,Y,TK_end)
    zlim([-0.2, 0.2]) 
    title('Temperature')
    subplot(2,5,10) % plot absolute error between exact and computed temperature at final time
    surf(X,Y,abs(TK_end-T_end)) %zlim([0, 3])
    title('Abs. temp. error')
    pause(0.01)
    
    % display
    if flag ~= 2
        disp(['' num2str(i) ') res = ' num2str(res(i)) ', DP = ' num2str(DP) ', rel = ' num2str(rel(i)) ', step(LMM) = ' num2str(step) ' '])
        %disp(['       (K11) err = ' num2str(ERR(i,1)) ', err_int = ' num2str(ERR(i,2)) ''])
        %disp(['       (K22) err = ' num2str(ERR(i,3)) ', err_int = ' num2str(ERR(i,4)) ''])
    end
    
    % update counter
    i = i + 1;
end

%% DISPLAYING RESULTS (RECONSTRUCTED CONDUCTIVITY)
K11 = Ks(1:(n+1)^2);
K11 = reshape(K11,n+1,n+1); K11 = K11(:);
K22 = Ks((n+1)^2+1:end);
K22 = reshape(K22,n+1,n+1); K22 = K22(:);
err1 = norm(Ke11(:) - K11)/norm(Ke11(:));
err2 = norm(Ke22(:) - K22)/norm(Ke22(:));
Ke11_int = Ke11(2:end-1,2:end-1);
Ke22_int = Ke22(2:end-1,2:end-1);
K11_full = reshape(K11,n+1,n+1);
K22_full = reshape(K22,n+1,n+1);
K11_int = K11_full(2:end-1,2:end-1);
K22_int = K22_full(2:end-1,2:end-1);
err_int1 = norm(Ke11_int(:) - K11_int(:))/norm(Ke11_int(:));
err_int2 = norm(Ke22_int(:) - K22_int(:))/norm(Ke22_int(:));

% display
disp(['(K11) error: ' num2str(err1) ', int. error: ' num2str(err_int1) ', iter: ' num2str(i-1) ''])
disp(['(K22) error: ' num2str(err2) ', int. error: ' num2str(err_int2) ''])
%disp('-------------------------------------------------------------')

sol(6:10) = [err1,err_int1,err2,err_int2,i-1]; % save data
if flag == 2 % DP stopping criteria
    sol(12) = 1;
elseif flag == 3 % relative error stopping
    sol(12) = 0;
else % kmax stopping
    sol(12) = 2;
end

% plot (approximate conductivities)
figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;
for jj = 1:2
    if jj == 1 % K11
        K = Ks(1:(n+1)^2);
        K = reshape(K,n+1,n+1); K = K(:);
        Ke = Ke11;
        ii = 0; % counter for plotting
    elseif jj == 2 % K22
        K = Ks((n+1)^2+1:end);
        K = reshape(K,n+1,n+1); K = K(:);
        Ke = Ke22;
        ii = 4; % counter for plotting
    end
    subplot(2,4,ii+2)
    contour(X,Y,reshape(Ke,n+1,n+1),20);
    axis square
    title('exact (contour)')
    subplot(2,4,ii+3)
    contour(X,Y,reshape(K,n+1,n+1),20);
    axis square
    %title(['residual = ' num2str(0.5*norm(UK_pos-Up)^2) ', DP = ' num2str(DP) ' '])
    title('approx. (contour)')
    subplot(2,4,ii+4)
    %surf(X,Y,reshape(Ks,n+1,n+1),'facecolor','white','edgecolor','interp'), zlim([0.05, 0.3]) 
    surf(X,Y,reshape(K,n+1,n+1),'facecolor','white','edgecolor','interp','linewidth',1.5), 
    zlim([0, 3]) 
    %zlim([0.08, 0.15])
    %zlim([0.08, 0.18])
    %zlim([0.05, 0.18])
    axis square
    %title(['error = ' num2str(error) ', int. error = ' num2str(err_int) ' '])
    title('approx.')
    subplot(2,4,ii+1)
    %surf(X,Y,reshape(Ks,n+1,n+1),'facecolor','white','edgecolor','interp'), zlim([0.05, 0.3]) 
    surf(X,Y,reshape(Ke,n+1,n+1),'facecolor','white','edgecolor','interp','linewidth',1.5), 
    zlim([0, 3])
    axis square
    title('exact')
end

%ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
%ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot (approximate conductivities)
figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;
for jj = 1:2
    if jj == 1 % K11
        K = Ks(1:(n+1)^2);
        K = reshape(K,n+1,n+1); K = K(:);
        Ke = Ke11;
        ii = 0; % counter for plotting
    elseif jj == 2 % K22
        K = Ks((n+1)^2+1:end);
        K = reshape(K,n+1,n+1); K = K(:);
        Ke = Ke22;
        ii = 4; % counter for plotting
    end
    subplot(2,4,ii+4)
    surf(X,Y,reshape(K,n+1,n+1),'facecolor','white','edgecolor','interp','linewidth',1.5), 
    zlim([0, 3]) 
    axis square
    %title('approx.')
end
ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);
% 
% export_fig cond_exact_n15_N10_delta001.eps -transparent

%% DISPLAYING RESULTS (TEMPERATURE RECONSTRUCTION)
T_end = reshape(Uex(:,end),n+1,n+1); % exact temperature at final time
TK = ort_generate_T(Ks,T0,mydata);
%TK = TK(:,2:end);
TK = TK(:,1 + (1:N)*(et+1)); % select time steps
TK_end = reshape(TK(:,end),n+1,n+1); % computed temperature at final time
Uex = Uex(:);
TK = TK(:);

errT = norm(TK-Uex)/norm(Uex); % error at all time steps
disp(['Temperature reconstruction error = ' num2str(errT) ''])
disp('-------------------------------------------------------------')

sol(11) = errT;

figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;

subplot(2,4,1) % plot computed temperature at final time
surf(X,Y,TK_end,'facecolor','white','edgecolor','interp','linewidth',1.5)
zlim([-0.2, 0.2]) 
%title(['Temperature at t = ' num2str(tt0(end)) ' '])
subplot(2,4,5) % plot absolute error between exact and computed temperature at final time
surf(X,Y,abs(TK_end-T_end),'facecolor','white','edgecolor','interp','linewidth',1.5), %zlim([0, 3])
%title('Absolute error')

ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);

%export_fig temp_exact_n15_N10_delta0001_2.eps -transparent












