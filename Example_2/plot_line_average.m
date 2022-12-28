% Code to plot lines of recovered and exact conductivity.

clear; close all;

% plot positioning
aa = (1920-1366)/1920*0.5; % notebook -> monitor display
bb = (1080-768)/1080*0.5;
%posi = [0 0 1 1]; % full screen (any monitor)
posi = [aa bb 1366/1920 768/1080]; % centered, notebook size

%% LOAD DATA
load('sol_f_k1_D2_norm2_over.mat')
n = 15;
repetition = size(Iter,2);
y1 = 0.2;
y2 = 0.8;
yline1 = y1*ones(n+1,1); % line position
yline2 = y2*ones(n+1,1);

[~,x] = cheby1(n,0,1);
y = x;
[X,Y] = ndgrid(x,y);
run('data.m') % original data
run('data_adequation.m') % data adequation to [0,1]x[0,1]
Ke11 = k11(X,Y); % exact K11
Ke22 = k22(X,Y); % exact K22
Ke11 = griddedInterpolant(X,Y,Ke11); % to guess values out of cheby's grid
Ke22 = griddedInterpolant(X,Y,Ke22);

%% LINE PLOTS
% plot for k11, in line, for NL 1.5% and different measurement scenarios

figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;

% plot exact
PlotOptions.LineWidth = 1.2;
PlotOptions.LineStyle = '-';
PlotOptions.Color = 'r';
PlotOptions.Marker = 'none';
subplot(2,4,1)
%plot(x,Ke11(x,yline1),PlotOptions), hold on
plot(x,Ke11(x,0.2*ones(n+1,1)),PlotOptions), hold on
subplot(2,4,2)
plot(x,Ke11(x,0.4*ones(n+1,1)),PlotOptions), hold on
%plot(x,Ke22(x,yline1),PlotOptions), hold on
subplot(2,4,3)
plot(x,Ke11(x,0.6*ones(n+1,1)),PlotOptions), hold on
subplot(2,4,4)
plot(x,Ke11(x,0.8*ones(n+1,1)),PlotOptions), hold on
%plot(x,Ke11(x,yline2),PlotOptions), hold on

for i = 1:3
    K = zeros(2*(n+1)^2,1);
    for j = 1:repetition
        K = K + Ksol(3*(i-1)+2,:,j)';
    end
    K = K/repetition; % average
    K11 = K(1:(n+1)^2);
    K22 = K((n+1)^2+1:end);
    K11 = reshape(K11,n+1,n+1);
    K22 = reshape(K22,n+1,n+1);
    K11 = griddedInterpolant(X,Y,K11); % to guess values out of cheby's grid
    K22 = griddedInterpolant(X,Y,K22);
    
    % plot options
    PlotOptions.LineWidth = 1.2;
    if i == 1 % NL = 0.02
        PlotOptions.LineStyle = '-';
        PlotOptions.Color = 'b';
        PlotOptions.Marker = 'o';
        PlotOptions.MarkerSize = 6;
    elseif i == 2 % NL = 0.04
        PlotOptions.LineStyle = '-';
        PlotOptions.Color = 'k';
        PlotOptions.Marker = '.';
        PlotOptions.MarkerSize = 14;
    elseif i == 3 % NL = 0.06
        PlotOptions.LineStyle = '-';
        PlotOptions.Color = [0 0.7 0.3];
        PlotOptions.Marker = 'v';
        PlotOptions.MarkerSize = 6;
    end

    % plots
    subplot(2,4,1)
    %plot(x,K11(x,yline1),PlotOptions), hold on
    plot(x,K11(x,0.2*ones(n+1,1)),PlotOptions), hold on
    subplot(2,4,2)
    %plot(x,K22(x,yline1),PlotOptions), hold on
    plot(x,K11(x,0.4*ones(n+1,1)),PlotOptions), hold on
    %plot(x,K22(x,yline1),PlotOptions), hold on
    subplot(2,4,3)
    %plot(x,K11(x,yline2),PlotOptions), hold on
    plot(x,K11(x,0.6*ones(n+1,1)),PlotOptions), hold on
    subplot(2,4,4)
    %plot(x,K11(x,yline2),PlotOptions), hold on
    plot(x,K11(x,0.8*ones(n+1,1)),PlotOptions), hold on
    
    
    %subplot(2,4,4)
    %plot(x,K22(x,yline2),PlotOptions), hold on
    %plot(x,K11(x,0.8*ones(n+1,1)),PlotOptions), hold on
    
end

% adjust plots
for i = 1:4
    subplot(2,4,i)
    axis square
    if i == 1
        legend('Exact','Full grid','MS1','MS2','location','northwest')
        xlabel('x')
        title(['k_{11}(x,' num2str(0.2) ')'])
        %ylim([0.1 0.24])
    elseif i == 2
        xlabel('x')
        title(['k_{11}(x,' num2str(0.4) ')'])
        %ylim([0.11 0.26])
    elseif i == 3
        xlabel('x')
        title(['k_{11}(x,' num2str(0.6) ')'])
        %ylim([0.14 0.28])
        %ylim([0.12 0.25])
        %yticks([0.12, 0.15, 0.2, 0.25])
    elseif i == 4
        xlabel('x')
        title(['k_{11}(x,' num2str(0.8) ')'])
        %ylim([0.14 0.26])
        %yticks([0.14:0.04:0.26])
    end
    %legend('Exact','NL=0.1%','NL=1.5%','NL=3%','location','northwest')
    ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
    ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);
end


























