% Routine to treat results for repetition tests (LMM and lsqnonlin), for
% the 2-norm tests (generally using relative residual criterion), including
% upper and lower bounds to lsqnonlin's solver.

% EXPERIMENTAL DATA CASE.

clear; close all;

%% PARAMETERS
chosenNorm = 2; % select norm
%final_time = 60; % final time
%time_steps = 180; % (N) number of time steps
K0 = 55; % initial guess
%delta = '0'; % noise level in the performed trials

% plot positioning
aa = (1920-1366)/1920*0.5; % notebook -> monitor display
bb = (1080-768)/1080*0.5;
%posi = [0 0 1 1]; % full screen (any monitor)
posi = [aa bb 1366/1920 768/1080]; % centered, notebook size

%% TREAT RESULTS
%storage = zeros(3,3); % storage results: [errorK, RTRE, max(iter)]

% load data
load(sprintf('EXPrep10k%dt60N180norm2delta0reg1.mat',K0))
    
% parameters
n = parameters(1); % number of Chebyshev's points
N = parameters(2); % number of time steps
x = X(:,1); % Chebyshev's grid (in [0,1], x direction)
rep = length(iter_LMM); % number of repetitions
ell = 12; % total number of sensors (spatial measurements)
l1 = 0.1;
l2 = 0.015;

% empty space to store results
iter = zeros(3,3); % iterations [min, max, mean]
nK = zeros(rep,3); % relative error for K
nPT = zeros(rep,3); % relative partial temperature error (RTRE)

for ii = 1:3
    if ii == 1 % LMM case
        K = K_LMM;
        TK = TK_LMM;
    elseif ii == 2 % TRR case
        K = K_TRR;
        TK = TK_TRR;
    elseif ii == 3 % TRR case (upper and lower bounds)
        K = K_TRR_ul;
        TK = TK_TRR_ul;
    end
    
    T_mean = zeros(ell,N);
    for i = 1:rep

        % conductivity errors
        nK(i,ii) = norm(K(:,i)-Ke,chosenNorm)/norm(Ke,chosenNorm);

        % temperature inner errors
        T = TK(:,i);
        T = reshape(T,(n+1)^2,N);
        qT = zeros(ell,N);
        for jj = 1:N
            qT(:,jj) = interp2(l1*X',l2*Y',reshape(T(:,jj),n+1,n+1)',xs',ys');
        end
        T = qT;
        nPT(i,ii) = norm(T(:)-U(:),chosenNorm)/norm(U(:),chosenNorm);
        
        T_mean = T_mean + T;

    end
    
    % save temperature mean values
    T_mean = T_mean/rep; % arithmetic mean
    if ii == 1
        T_LMM = T_mean;
    elseif ii == 2
        T_TRR = T_mean;
    elseif ii == 3
        T_TRR_ul = T_mean;
    end

end

iterations = [iter_LMM, iter_TRR, iter_TRR_ul]; % iterations

% results (for LMM (first line) and TRR (second line))
results = zeros(3,3);
results(:,1) = mean(nK)'; % mean of conductivity error
results(:,2) = mean(nPT)'; % mean of partial temperature error
results(:,3) = max(iterations)'; % max iterations

% store results
storage = results;

% build table
T = array2table(storage,'VariableNames',{'errK','RTRE','MI'},'RowNames',{'LMM',...
    'TRR','TRR ul'}); % table of results

% display results
disp(' ')
disp('Matrix with results:')
disp(' ')
disp(storage)
disp(' ')
disp(T)


%% PLOTS (another selection)
figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;

% plot position
subplot(2,3,1)

% plot
x = X(:,1); % Chebyshev's grid (in [0,1], x direction)
plot(x,Ke,'r','linewidth',1.3), hold on,
plot(x,mean(K_LMM,2),'b-o','linewidth',1.3),
plot(x,mean(K_TRR_ul,2),'k.-','linewidth',1.3,'MarkerSize',10,'color',[0 0 0])
%plot(x,mean(K_TRR_ul,2),'-v','linewidth',1.3,'color',[0 0.5 0])
%xlim([1 n+1])
%ylim([0 100])
ylim([0 100])
legend('Exact','LMM','lsqnonlin','location','northwest')
xlabel('x')
ylabel('Cond. (W/m°C)')
axis square

ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);

%% PLOT TEMPERATURE ERRORS AT SELECTED THERMOCOUPLES
figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;

tt0 = t0(2:end);

thermocouples = [1, 7, 8, 12]; % chosen thermocouples for plotting
places = [1, 2, 3, 4]; % positions for subplots

for i = 1:4
    
    ii = thermocouples(i); % chosen thermocouple
    
    subplot(3,2,places(i))
    %subplot('position',[0.1 0.5 0.65 0.2])
    selected = [1:3:180, 180];
    
    % plot
    %plot(tq0,Tq(ii,:),'r','linewidth',1.3), hold on,
    %plot(tt0,U(ii,:),'r','linewidth',1.3), hold on,
    %plot(tt0,T_TRR(ii,:),'-','linewidth',1.3,'color',[0 0.5 0]),
    %plot(tt0,T_LMM(ii,:),'-','linewidth',1.3,'color',[0 0 1]),
    plot(tt0(selected),abs(T_LMM(ii,selected)-U(ii,selected))./U(ii,selected),'-','linewidth',1.3,'color',[0 0 1]), hold on,
    %plot(tt0,T_TRR_ul(ii,:),'-','linewidth',1.3,'MarkerSize',10,'color',[0 0 0]),
    plot(tt0(selected),abs(T_TRR_ul(ii,selected)-U(ii,selected))./U(ii,selected),'--','linewidth',1.3,'MarkerSize',10,'color',[0 0 0])
    
    ylim([0 0.15])
    xlabel('Time (s)')
    ylabel('Relative error')
    title(sprintf('Thermocouple T%d',ii))
    legend('LMM','lsqnonlin','location','northwest')

end

ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);










































