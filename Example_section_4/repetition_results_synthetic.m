% Routine to treat results for repetition tests (LMM and lsqnonlin), for
% the 2-norm tests (generally using relative residual criterion), including
% upper and lower bounds to lsqnonlin's solver.

% SYNTHETIC DATA CASE.

clear; close all;

%% PARAMETERS
chosenNorm = 2; % select norm
K01 = 25; % first initial guess
K02 = 55; % second initial guess
delta = '0015'; % noise level in the performed trials

% plot positioning
aa = (1920-1366)/1920*0.5; % notebook -> monitor display
bb = (1080-768)/1080*0.5;
%posi = [0 0 1 1]; % full screen (any monitor)
posi = [aa bb 1366/1920 768/1080]; % centered, notebook size

%% TREAT RESULTS
storage = zeros(18,3); % storage results: [errorK, RTRE, max(iter)]

for situation = 1:6
    
    % t60N60
    if situation == 1
        load(sprintf('ulrep10k%dt60N60norm2delta%sreg1.mat',K01,delta))
    elseif situation == 2
        load(sprintf('ulrep10k%dt60N60norm2delta%sreg1.mat',K02,delta))
    
    % t30N30
    elseif situation == 3
        load(sprintf('ulrep10k%dt30N30norm2delta%sreg1.mat',K01,delta))
    elseif situation == 4
        load(sprintf('ulrep10k%dt30N30norm2delta%sreg1.mat',K02,delta))
        
    %t30N60
    elseif situation == 5
        load(sprintf('ulrep10k%dt30N60norm2delta%sreg1.mat',K01,delta))
    elseif situation == 6
        load(sprintf('ulrep10k%dt30N60norm2delta%sreg1.mat',K02,delta))
    end
    
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

        end
        
    end

    iterations = [iter_LMM, iter_TRR, iter_TRR_ul]; % iterations

    % results (for LMM (first line) and TRR (second line))
    results = zeros(3,3);
    results(:,1) = mean(nK)'; % mean of conductivity error
    results(:,2) = mean(nPT)'; % mean of partial temperature error
    results(:,3) = max(iterations)'; % max iterations
    
    % store results
    storage(3*(situation-1)+1:3*situation,:) = results;

end

% build table
storage = [storage(1:6,:), storage(7:12,:), storage(13:18,:)];
T = array2table(storage,'VariableNames',{'errK_0','RTRE_0','MI_0','errK_1',...
    'RTRE_1','MI_1','errK_2','RTRE_2','MI_2'},'RowNames',{'(k01) LMM',...
    '(k01) TRR','(k01) TRR ul','(k02) LMM','(k02) TRR','(k02) TRR ul'}); % table of results

% display results
disp(' ')
disp('Matrix with results:')
disp(' ')
disp(storage)
disp(' ')
disp('------------------------------------')
disp('Legend for the table:')
disp('errK_0: case t60, N60')
disp('errK_1: case t30, N30')
disp('errK_2: case t30, N60')
disp('------------------------------------')
disp(' ')
disp(T)

%% PLOTS
    
for situation = 1:3
    
    if situation == 1 % t60N60
        load(sprintf('ulrep10k%dt60N60norm2delta%sreg1.mat',K02,delta))
    elseif situation == 2 % t30N30
        load(sprintf('ulrep10k%dt30N30norm2delta%sreg1.mat',K02,delta))
    elseif situation == 3 % t30N60
        load(sprintf('ulrep10k%dt30N60norm2delta%sreg1.mat',K02,delta))
    end
    
    % plot position
    figure,
    fig = gcf;
    fig.Units = 'normalized';
    fig.OuterPosition = posi;
    subplot(2,3,situation)
    
    % plot
    x = X(:,1); % Chebyshev's grid (in [0,1], x direction)
    plot(x,Ke,'r','linewidth',1.3), hold on,
    plot(x,mean(K_LMM,2),'b-o','linewidth',1.3),
    plot(x,mean(K_TRR_ul,2),'k.-','linewidth',1.3,'MarkerSize',10,'color',[0 0 0])
    %plot(x,mean(K_TRR_ul,2),'-v','linewidth',1.3,'color',[0 0.5 0])
    %xlim([1 n+1])
    %ylim([0 100])
    ylim([0 100])
    if situation == 1
        %title('k_0 \approx 2')
%         legend('Exact','lsqnonlin (DP)','lsqnonlin (Rel.)','LMM (DP)','LMM (Rel.)','location','best')
        legend('Exact','LMM','lsqnonlin','location','northwest')
    end
    xlabel('x')
    ylabel('Cond. (W/m°C)')
    axis square
    
    ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
    ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);

end


















