% Code to analyze results calculated in repetition tests (essentially, 
% change the norm chosed, but other computations are certainly allowed as 
% well, if necessary).

clear; close all;

%% Parameters
%load('sol_selection_over13.mat') % load repetition data
load('sol_cao_1over4_D1_norm2_over_v2.mat')
% sol = [points,delta,err1,err_int1,err2,err_int2,iter_min,iter_max,errT,errTinner,success]

n = 15; % number of Cheby's points
N = 10; % time steps
tf = 1; % final time
repetition = 30; % number of repetitions
chosenNorm = 2; % chosen norm to be used

t0 = linspace(0,1,N+1);
[~,x] = cheby1(n,0,1); % Cheby's points
y = x;
[X,Y] = ndgrid(x,y); % grid
run('data.m') % load problem's data
run('data_adequation.m')

% exact conductivity
Ke11 = k11(X,Y);
Ke11 = Ke11(:);
Ke22 = k22(X,Y);
Ke22 = Ke22(:);
Ke = [Ke11(:); Ke22(:)]; 

% exact temperature
Uex = zeros((n+1)^2,N);
for i = 1:N
    TT = T(X,Y,t0(i+1)); % calculate exact temperature at t_{i+1}
    Uex(:,i) = TT(:); % store exact temperature
end
Te = Uex(:); % vector version

%% Compute quantities

S = zeros(size(sol,1), 11);  % [points,delta,errK,err1,err2,iterMin,iterMax,iterMean,errT,errorPT,success]
S(:,[1,2,6,7,11]) = sol(:,[1,2,7,8,11]); % allocate some of the data
for i = 1:9
        
    % determine inner points
    if sol(i,1) == 1  % 'full'
        ell_inner = 1:(n+1)^2;
        ell_inner = ell_inner';
    elseif sol(i,1) == 2  % 'inner'
        xi = zeros(n+1,1);
        yi = zeros(n+1,1);
        for ii = 1:n+1
            if l1*x(ii) >= x1 && l1*x(ii) <= x2
                xi(ii) = 1;
            end
            if l2*y(ii) >= y1 && l2*y(ii) <= y2
                yi(ii) = 1;
            end
        end
        ell_inner = zeros((n+1)^2,1); % to select inner Jacobian lines
        for ii = 1:n+1
            for j = 1:n+1
                if xi(j) == 1 && yi(ii) == 1 
                    ell = (j-1) + (ii-1)*(n+1) + 1;
                    ell_inner(ell) = ell;
                end
            end
        end
        ell_inner = nonzeros(ell_inner);
    elseif sol(i,1) == 3  % 'selection'
        ell_inner = [51,54,57,60,84,87,90,93,117,120,123,126,147,150,153,156,180,183,186,189]';
    end
    ell_inner = sort(ell_inner); % to reorder the points
    
    %Ks = Ksol(i,:,:); % conductivities
    %TKs = TKsol(i,:,:); % temperature values
    Tep = Uex(ell_inner,:); % inner points exact temperature
    Tep = Tep(:);
    results = zeros(repetition, 6); % store results
    for k = 1:repetition
        
        % conductivity
        K = Ksol(i,:,k);
        K11 = K(1:(n+1)^2);
        K22 = K((n+1)^2+1:end);
        errK = norm(K(:)-Ke, chosenNorm)/norm(Ke, chosenNorm);
        err1 = norm(K11(:)-Ke11, chosenNorm)/norm(Ke11, chosenNorm);
        err2 = norm(K22(:)-Ke22, chosenNorm)/norm(Ke22, chosenNorm);

        % temperature
        TK = TKsol(i,:,k);
        errT = norm(TK(:)-Te, chosenNorm)/norm(Te, chosenNorm);
        TK = reshape(TK,(n+1)^2,N);
        TKp = TK(ell_inner,:);
        TKp = TKp(:);
        errPT = norm(TKp-Tep, chosenNorm)/norm(Tep, chosenNorm);
        
        % iteration
        iter = Iter(i,k,5);
        
        % store
        results(k,:) = [errK, err1, err2, iter, errT, errPT];
        
    end
    
    % compute means
    S(i,[3,4,5,8,9,10]) = mean(results);
    
end

S(:,:) = S([1 4 7 2 5 8 3 6 9],:); % reorder to noise levels

disp('Solution: ')
disp(S)



























