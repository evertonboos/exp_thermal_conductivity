% Implementation and solution using pseudospectral method for
% discretization and Crank-Nicolson for the resulting system of ODE,
% related to the heat flow equation. Here, we aim to solve the inverse
% problem of finding K given the 'exact' temperature and an approximation
% of the initial K, by Levenberg-Marquardt's method.

% EXPERIMENTAL DATA CASE, SYNTHETIC EXPERIMENTS.

%clear; close all;

%% PARAMETERS AND GRID DETAILS
n = 15; % Chebyshev's points 
N = 60; % (30,60) number of equally spaced intervals on time discretization
tf = 60; % (30, 60) final time
K0 = 55; % initial value of k
global CRITERION, CRITERION = 1; % '0' = DP, '1' = relative residual
global chosenNorm, chosenNorm = 2; % '2' for euclidian norm, 'inf' for infinity norm
delta = 0.015; % noise level
repetition_number = 1; % number of tests repetition; used '10' in the paper
lowerBound = 20; % lower bound for lsqnonlin
upperBound = 60; % upper bound for lsqnonlin

% saveFileName = sprintf('ulrep10k%dt%dN%dnorm2delta0015reg1.mat',K0,final_t,N); % string to be used as name for the .mat file containing the results of this code

tau = 1.1; % DP parameter
kmax = 50; % maximum number of iterations
alpha = 0; % influence of conductivity on the boundary

parameters = [n,N,tf,delta,lowerBound,upperBound]; % save data

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
%inds = [1:8, 12];
inds = 1:12;
xs = xs(inds);
ys = ys(inds);
ell = length(xs);

% exact conductivity
Ke = k(x,y);
Ke = Ke(:);

%% TEMPERATURE MEASUREMENTS AND NOISE
load('sensor_n40_N600.mat') % measurements on sensor locations
nd = size(Tq,2); % number of measured instances
%delta_t = tf/nd; % time delta
tt0 = zeros(1,N+1);

% time steps determination:
delta_nd = floor( (tf/tq0(end)) * nd/N ); % to choose time steps 'equally spaced'
%delta_nd = 1;

U = zeros(ell,N); % store experimental data (observe initial condition is not included)

for i = 1:N
    instance = i*delta_nd;
    tt0(i+1) = tq0(instance); % time instances
    TT = Tq(inds,instance); % selection of sensor values at correct time instance
    U(:,i) = TT(:);
end

%tt0 = linspace(0,tf,N+1);

t0 = tt0;
mydata = struct('cheby', {n,N,t0,tt0,et,alpha,beta,D,x,y,P,hat_D,d0,dn,d0_e1,dn_en,X,Y});

T_end = U(:,end); % final temperature

rng('default'); % to set tests' same condition always
K_initial = K0 + randn(numel(Ke),repetition_number); % initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Noise = randn(ell*N,repetition_number); % noise computation

% store values variables
K_LMM = zeros(numel(Ke),repetition_number); % solutions (LMM)
TK_LMM = zeros(N*(n+1)^2,repetition_number); % temperature (LMM)
iter_LMM = zeros(repetition_number,1); % iterations (LMM)
K_TRR = zeros(numel(Ke),repetition_number); % solutions (trust-region-reflective)
TK_TRR = zeros(N*(n+1)^2,repetition_number); % temperature (trust-region-reflective)
iter_TRR = zeros(repetition_number,1); % iterations (trust-region-reflective)
K_TRR_ul = zeros(numel(Ke),repetition_number); % solutions (trust-region-reflective WITH UPPER AND LOWER BOUNDS)
TK_TRR_ul = zeros(N*(n+1)^2,repetition_number); % temperature (trust-region-reflective WITH UPPER AND LOWER BOUNDS)
iter_TRR_ul = zeros(repetition_number,1); % iterations (trust-region-reflective WITH UPPER AND LOWER BOUNDS)

for ii = 1:repetition_number
    
    % display repetition
    disp(['REPETITION: ' num2str(ii) ''])
    
    % noise introduction
    nU = norm(U(:),chosenNorm);
    %noise = randn(ell*N,1);
    noise = Noise(:,ii);
    noise = noise/norm(noise,chosenNorm);
    errU = delta*nU; 
    Up = U + delta*nU*reshape(noise,ell,N); % perturbed U, ||Up-U|| = delta*||U|| (as vector)
    Up = Up(:);

    % stopping criteria setting
    global DP
    DP = tau*errU; % DP criterion

    %% CPS + LEVENBERG-MARQUARDT'S METHOD

    gamma = 0.7; % LMM constant, \gamma \in (0,1)
    sigma = 1e-2; % Armijo's constant, \sigma \in (0,1)

    % first iterate K
    Ks = K_initial(:,ii);

    % LMM scaling matrix
    L = get_l(n+1,1); % first or second order discrete derivative
    LL = [kron(eye(n+1),L); kron(L,eye(n+1))];
    D = L; % scaling with first or second order derivative (isox)

    % storage errors and residuals
    res = zeros(kmax,1); 
    rel = res;
    Jcond = res; % conditioning of [J,L]
    ERR = zeros(kmax,4); % [err(K11), int_err(K11), err(K22), int_err(K22)]
    % storage errors and residuals

    % ITERATIONS
    Kc = zeros(numel(Ke),kmax+1); % storage computed solutions
    Kc(:,1) = Ks;
    flag = 0;
    i = 1; % counter
    while (i <= kmax && flag == 0)
        JJ = isox_Jacobian(Ks,mydata); % Jacobian at time steps t0
        J = sparse(ell*N,numel(Ke)); % inner Jacobian
        for j = 1:N-1
            Jk = zeros(ell,numel(Ke)); % space for Jacobian at every time step
            for jj = 1:numel(Ke)
                Jk(:,jj) = interp2(l1*X',l2*Y',reshape(JJ(:,jj,j),n+1,n+1)',xs',ys'); % interpolate Jacobian relative to sensor locations
            end
            int = (j-1)*ell+1:j*ell;
            J(int,:) = Jk; 
        end

        TK_ant = isox_generate_T(Ks,T0,mydata); % solution for K
        %TK_ant = TK_ant(:,2:end);
        TK_ant = TK_ant(:,1 + (1:N)*(et+1)); % select time steps
        qTK = zeros(ell,N);
        for jj = 1:N
            qTK(:,jj) = interp2(l1*X',l2*Y',reshape(TK_ant(:,jj),n+1,n+1)',xs',ys');
        end
        TK_ant = qTK;
        %TK_ant = get_inner(TK_ant,n,xi,yi); % get inner points
        UK_ant = TK_ant(:); % vectorized solution for K
        F_ant = norm(UK_ant-Up);
        mu = F_ant^2;
        
        % compute direction d_k (LMM)
        d = ( J'*J + mu*(D'*D) ) \ ( -J'*(UK_ant-Up) ); % direction d_k

        TK_pos = isox_generate_T(Ks+d,T0,mydata);
        TK_end = TK_pos(:,end);
        %TK_pos = TK_pos(:,2:end);
        TK_pos = TK_pos(:,1 + (1:N)*(et+1)); % select time steps
        qTK = zeros(ell,N);
        for jj = 1:N
            qTK(:,jj) = interp2(l1*X',l2*Y',reshape(TK_pos(:,jj),n+1,n+1)',xs',ys');
        end
        TK_pos = qTK;
        %TK_pos = get_inner(TK_pos,n,xi,yi); % get inner points
        UK_pos = TK_pos(:);
        F_pos = norm(UK_pos-Up);

        % New Iterate Definition (LMM with line search)
        if F_pos <= gamma*F_ant % new iterate (without line search)
            Ks = Ks + d;
            step = 1; 
        elseif flag == 0 % LINE SEARCH
            m = 1;
            TK = isox_generate_T(Ks + gamma^m*d,T0,mydata);
            %TK = TK(:,2:end);
            TK = TK(:,1 + (1:N)*(et+1)); % select time steps
            qTK = zeros(ell,N);
            for jj = 1:N
                qTK(:,jj) = interp2(l1*X',l2*Y',reshape(TK(:,jj),n+1,n+1)',xs',ys');
            end
            TK = qTK;
            %TK = get_inner(TK,n,xi,yi); % get inner points
            UK = TK(:);
            F = norm(UK-Up);
            gradphi_d = (J'*(UK_ant-Up))'*d; 
            while 0.5*F^2 - 0.5*F_ant^2 > sigma*gamma^m*gradphi_d % Armijo
                if gamma^m < 1e-10
                    disp('Armijo failed, aborting...')
                    flag = 1; % abortion flag
                    break
                end
                m = m+1;
                TK = isox_generate_T(Ks + gamma^m*d,T0,mydata);
                TK_end = TK(:,end);
                %TK = TK(:,2:end);
                TK = TK(:,1 + (1:N)*(et+1)); % select time steps
                qTK = zeros(ell,N);
                for jj = 1:N
                    qTK(:,jj) = interp2(l1*X',l2*Y',reshape(TK(:,jj),n+1,n+1)',xs',ys');
                end
                TK = qTK;
                %TK = get_inner(TK,n,xi,yi); % get inner points
                UK = TK(:);
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

        % Stopping criteria

        res(i) = norm(UK_pos-Up,chosenNorm); % residual

        % relative residual test
        if CRITERION == 1
            if i >= 2
                rel(i) = abs(res(i) - res(i-1))/res(i-1);
                %rel(i) = norm(Kc(:,i)-Kc(:,i-1))/norm(Kc(:,i));
                if rel(i) <= 1e-3
                    disp('-------------------------------------------------------------')
                    disp(['Relative residual criteria reached! Iterations: ' num2str(i) ' '])
                    flag = 3;
                end
            end
        end
        
        % DP test
        if CRITERION == 0
            if res(i) <= DP
                disp('-------------------------------------------------------------')
                disp(['DP reached! Iterations: ' num2str(i) ' '])
                flag = 2;
            end
        end

        disp(['' num2str(i) ') res = ' num2str(res(i)) ', DP = ' num2str(DP) ', rel = ' num2str(rel(i)) ' '])

        % update counter
        i = i + 1;
    end
    
    % Store
    K_LMM(:,ii) = Ks; % solution
    iter_LMM(ii) = i-1; % iterations
    TK = isox_generate_T(Ks,T0,mydata);
    TK = TK(:,2:end);
    TK_LMM(:,ii) = TK(:); % temperature

    %% CPS + lsqnonlin (trust-region-reflective)
    
    global sK sRes sNorm sGnorm % variables to store values from lsqnonlin
    sK = zeros(n+1,kmax+1);
    sRes = zeros(kmax+1,1);
    sNorm = zeros(kmax+1,1);
    sGnorm = zeros(kmax+1,1);

    k0 = K_initial(:,ii); % initial guess for minimization
    sK(:,1) = k0;

    % Inverse problem solution
    OPTIONS = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective',...
        'Display','off','FunctionTolerance',1e-6,'OutputFcn',@outDP,...  % 'Display','iter-detailed'
        'MaxIterations',kmax); %,'PlotFcn','optimplotx');

    [K,~,~,~,~,~,~] = ...
         lsqnonlin(@(K) iso_resi(K,t0,Up,ell,l1,l2,X,Y,xs,ys),k0,0*ones(size(k0)),[],OPTIONS); 
         %lsqnonlin(@(K) iso_resi(K,t0,Up,ell,l1,l2,X,Y,xs,ys),k0,0*ones(size(k0)),90*ones(size(k0)),OPTIONS); % TRR with lower and upper bounds
         %iso_resi(K,t0,Upp2,ell_inner,x,xs,ys),k0,zeros(size(k0)),[],OPTIONS);
         %lsqnonlin(@(K) iso_resi(K,t0,Upp2,ell_inner,x,xs,ys),k0,OPTIONS);

    % Store
    K_TRR(:,ii) = K; % solution
    nn = nnz(sRes); % iterations
    iter_TRR(ii) = nn;
    TK = iso_gen_T(K,t0); % temperature
    TK = TK(:,2:end);
    TK_TRR(:,ii) = TK(:);
    
    %% CPS + lsqnonlin (trust-region-reflective WITH UPPER AND LOWER BOUNDS)
    
    % global sK sRes sNorm sGnorm % variables to store values from lsqnonlin
    sK = zeros(n+1,kmax+1);
    sRes = zeros(kmax+1,1);
    sNorm = zeros(kmax+1,1);
    sGnorm = zeros(kmax+1,1);

    k0 = K_initial(:,ii); % initial guess for minimization
    sK(:,1) = k0;

    % Inverse problem solution
    OPTIONS = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective',...
        'Display','off','FunctionTolerance',1e-6,'OutputFcn',@outDP,...  % 'Display','iter-detailed'
        'MaxIterations',kmax); %,'PlotFcn','optimplotx');

    [K,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = ...
         lsqnonlin(@(K) iso_resi(K,t0,Up,ell,l1,l2,X,Y,xs,ys),k0,lowerBound*ones(size(k0)),upperBound*ones(size(k0)),OPTIONS); % TRR with lower and upper bounds

    % Store
    K_TRR_ul(:,ii) = K; % solution
    nn = nnz(sRes); % iterations
    iter_TRR_ul(ii) = nn;
    TK = iso_gen_T(K,t0); % temperature
    TK = TK(:,2:end);
    TK_TRR_ul(:,ii) = TK(:);
    
end

% % SAVE DATA
% save(saveFileName,'K_LMM','K_TRR','K_TRR_ul',...
%     'K_initial','Ke','TK_LMM','TK_TRR','TK_TRR_ul','U','Up','X','Y',...
%     'iter_LMM','iter_TRR','iter_TRR_ul','parameters','t0','xs','ys');











% End of main code

%% Auxiliar functions (for lsqnonlin)
% output function to introduce DP within iterations
function stop = outDP(x,optimValues,state)
stop = false;
global CRITERION chosenNorm % variable to select stopping criterion
global DP % load global variable 'DP', since it is not introduced in the function's parameters
global sK sRes sNorm sGnorm % to store values

ITER = optimValues.iteration + 1;
residual = optimValues.residual;
sK(:,ITER) = x;
sRes(ITER) = norm(residual,chosenNorm);
sNorm(ITER) = norm(x);
sGnorm(ITER) = norm(optimValues.gradient);

%disp('residual norm: ')
%disp(sqrt(optimValues.resnorm))

% Check if DP is satisfied
% 'resnorm' is the residual norm, i.e resnorm = f(K) = ||T(K)-Upp||^2
if CRITERION == 0
    if sRes(ITER) <= DP % test: ||T(K)-Upp|| < DP
        stop = true;
        disp('Stopped using DP')
        disp(ITER)
    end
end

%sgnorm = 0;
%snorm = 0;
sres = 0;

% Check for relative residual test
if CRITERION == 1
    if ITER >= 2
        sres = abs(sRes(ITER)-sRes(ITER-1))/abs(sRes(ITER-1)); % residual

        % Check if relative residual variation is small
        if (sres <= 1e-3 && sres ~= 0)
            stop = true;
            disp('Stopped using relative residual')
            disp(ITER)
        end

    end
end

disp(['' num2str(ITER) ') res = ' num2str(sRes(ITER)) ', DP = ' num2str(DP) ', rel = ' num2str(sres) ' '])

%disp([' ' num2str(ITER) ') residual = ' num2str(sRes(ITER)) ', Gnorm = ' num2str(sgnorm) ', ResNorm = ' num2str(sres) ', Knorm = ' num2str(snorm) ' '])

end
   
% Auxiliar function
function residual = iso_resi(K1,t0,Up,ell,l1,l2,X,Y,xs,ys)
% This function computes the solution to the direct problem using Approach
% 1 via conductivity K and time steps t0. Then, it
% computes the resitual norm || T(K)-Up ||, where Up corresponds to
% the input data (temperature measurements) provided in the whole Cheby's
% mesh. 'Up' can be given in matrix or vector form.

% Parameters
n = numel(K1) - 1;
N = length(t0) - 1;

% Direct problem solution
TK = iso_gen_T(K1,t0);
TK = TK(:,2:end); % remove initial condition
TKI = zeros(ell,N);
for jj = 1:N
    TKI(:,jj) = interp2(l1*X',l2*Y',reshape(TK(:,jj),n+1,n+1)',xs',ys');
end
%TKI = TK(ell_inner,:); % TO MATCH THIS TEST WITH THE APPROACH 1 TRIAL
%RES=(TKI(ell_inner,5:2:end)-Up(ell_inner,5:2:end));
%RES=(TKI(ell_inner,:)-Up(ell_inner,:));
RES=(TKI(:)-Up); % Up comparison on points of interest

% Residual norm
residual = RES(:);
end


% Auxiliar function
function TK = iso_gen_T(K1,t0)
% calculate temperature T at time steps given in t0 via K as vector or matrix
% containing the heat conductivity parameter (Approach 1).

% Parameters
n = numel(K1)-1; % determine 'n' from 'K'
K = K1*ones(1,n+1);
% Grid details
[D,x] = cheby1(n,0,1); % differentiation matrix 
y = x;
[~,P] = permuta(n); % permutation matrix
D1 = D(:,2:end-1);
D2 = D(2:end-1,:);
d0 = D(:,1); 
dn = D(:,n+1);
d0_e1 = sparse(n+1,n+1); 
d0_e1(:,1) = d0; % product d_0*e_{1}^T
dn_en = sparse(n+1,n+1); 
dn_en(:,end) = dn; % product d_n*e_{n+1}^T

[X,Y] = ndgrid(x,y);

% data from the problem
run('data.m')
run('data_adequation.m')

% Assembling coefficient matrix M
F = sparse((n+1)^2,(n+1)^2);
G = sparse((n+1)^2,(n+1)^2);

for j = 1:n+1
    int = (j-1)*(n+1)+1:j*(n+1);
    A = l1*h1(y(j))*d0_e1 + D1*diag(K(2:end-1,j))*D2 - l1*h2(y(j))*dn_en;
    B = l2*h3(x(j))*d0_e1 + D1*diag(K(j,2:end-1))*D2 - l2*h4(x(j))*dn_en;
    F(int,int) = [A(1,:); D2*diag(K(:,j))*D; A(end,:)];
    G(int,int) = [B(1,:); D2*diag(K(j,:))*D; B(end,:)];
end

M = (1/l1^2)*F + (1/l2^2)*P*G*P';

% IVP
N = length(t0)-1;
TK = zeros((n+1)^2,N+1);
TK(:,1) = T0; % initial condition

C = c(X,Y);
C = C(:);
C = sparse(1:(n+1)^2, 1:(n+1)^2, C); % 'C' matrix, diagonal sparse

for i = 1:N % Crank-Nicolson's method
    
    % t_{i+1}
    h = t0(i+1)-t0(i); % \Delta t
    
    Mm = C - (h/2)*M;
    Mp = C + (h/2)*M;
    
    % right-hand side (vector S(t))
    ua = zeros((n+1)^2,1); va = zeros((n+1)^2,1); % t_i (anterior)
    up = zeros((n+1)^2,1); vp = zeros((n+1)^2,1); % t_{i+1} (posterior)
    for j = 1:n+1
        int = (j-1)*(n+1)+1:j*(n+1);
        % t_i
        b = l1*(-d0*h1(y(j))*f1(y(j),t0(i)) + dn*h2(y(j))*f2(y(j),t0(i)));
        d = l2*(-d0*h3(x(j))*f3(x(j),t0(i)) + dn*h4(x(j))*f4(x(j),t0(i)));
        ua(int) = [b(1); zeros(n-1,1); b(end)];
        va(int) = [d(1); zeros(n-1,1); d(end)];
        % t_{i+1}
        b = l1*(-d0*h1(y(j))*f1(y(j),t0(i+1)) + dn*h2(y(j))*f2(y(j),t0(i+1)));
        d = l2*(-d0*h3(x(j))*f3(x(j),t0(i+1)) + dn*h4(x(j))*f4(x(j),t0(i+1)));
        up(int) = [b(1); zeros(n-1,1); b(end)];
        vp(int) = [d(1); zeros(n-1,1); d(end)];
    end
    Sa = (1/l1^2)*ua + (1/l2^2)*P*va + g(X,Y,t0(i)); % S(t_i) (anterior)
    Sp = (1/l1^2)*up + (1/l2^2)*P*vp + g(X,Y,t0(i+1)); % S(t_{i+1}) (anterior)
    
    % Crank-Nicolson's method
    TK(:,i+1) = Mm\( Mp*TK(:,i) + (h/2)*(Sa+Sp) );
end


end


  
  










