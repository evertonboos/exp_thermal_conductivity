% code to plot lines of conductivity.

clear; close all;

% plot positioning
aa = (1920-1366)/1920*0.5; % notebook -> monitor display
bb = (1080-768)/1080*0.5;
%posi = [0 0 1 1]; % full screen (any monitor)
posi = [aa bb 1366/1920 768/1080]; % centered, notebook size

%% LINE PLOT
%load('sol_cao_1over4.mat')
load('sol_cao_1over4_D1_norm2_over_v2.mat');
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

xx = X(:);
yy = Y(:);

% % K11 PLOTS
% figure,
% fig = gcf;
% fig.Units = 'normalized';
% fig.OuterPosition = posi;
% 
% % plot exact
% PlotOptions.LineWidth = 1.2;
% PlotOptions.LineStyle = '-';
% PlotOptions.Color = 'r';
% subplot(2,3,1)
% %plot(x,Ke11(x,yline1),PlotOptions), hold on
% plot(x,Ke11(x,0.2*ones(n+1,1)),PlotOptions), hold on
% subplot(2,3,2)
% plot(x,Ke11(x,0.5*ones(n+1,1)),PlotOptions), hold on
% %plot(x,Ke22(x,yline1),PlotOptions), hold on
% subplot(2,3,3)
% plot(x,Ke11(x,0.8*ones(n+1,1)),PlotOptions), hold on
% %plot(x,Ke11(x,yline2),PlotOptions), hold on
% 
% 
% %subplot(2,4,4)
% %plot(x,Ke11(x,0.8*ones(n+1,1)),PlotOptions), hold on
% %plot(x,Ke22(x,yline2),PlotOptions), hold on
% 
% for i = 1:3
%     K = zeros(2*(n+1)^2,1);
%     for j = 1:repetition
%         K = K + Ksol(i,:,j)';
%     end
%     K = K/repetition; % average
%     K11 = K(1:(n+1)^2);
%     K22 = K((n+1)^2+1:end);
%     K11 = reshape(K11,n+1,n+1);
%     K22 = reshape(K22,n+1,n+1);
%     K11 = griddedInterpolant(X,Y,K11); % to guess values out of cheby's grid
%     K22 = griddedInterpolant(X,Y,K22);
%     
%     % plot options
%     PlotOptions.LineWidth = 1.2;
%     if i == 1 % NL = 0.02
%         PlotOptions.LineStyle = '-';
%         PlotOptions.Color = 'k';
%         PlotOptions.Marker = '.';
%         PlotOptions.MarkerSize = 10;
%     elseif i == 2 % NL = 0.04
%         PlotOptions.LineStyle = '-';
%         PlotOptions.Color = 'b';
%         PlotOptions.Marker = 'o';
%         PlotOptions.MarkerSize = 6;
%     elseif i == 3 % NL = 0.06
%         PlotOptions.LineStyle = '-';
%         PlotOptions.Color = [0 0.8 0.5];
%         PlotOptions.Marker = '*';
%         PlotOptions.MarkerSize = 6;
%     end
% 
%     % plots
%     subplot(2,3,1)
%     %plot(x,K11(x,yline1),PlotOptions), hold on
%     plot(x,K11(x,0.2*ones(n+1,1)),PlotOptions), hold on
%     subplot(2,3,2)
%     %plot(x,K22(x,yline1),PlotOptions), hold on
%     plot(x,K11(x,0.5*ones(n+1,1)),PlotOptions), hold on
%     %plot(x,K22(x,yline1),PlotOptions), hold on
%     subplot(2,3,3)
%     %plot(x,K11(x,yline2),PlotOptions), hold on
%     plot(x,K11(x,0.8*ones(n+1,1)),PlotOptions), hold on
%     
%     
%     %subplot(2,4,4)
%     %plot(x,K22(x,yline2),PlotOptions), hold on
%     %plot(x,K11(x,0.8*ones(n+1,1)),PlotOptions), hold on
%     
% end
% 
% % adjust plots
% for i = 1:3
%     subplot(2,3,i)
%     axis square
%     if i == 1
%         legend('Exact','NL = 2%','NL = 4%','NL = 6%','location','northwest')
%         xlabel('x')
%         ylabel(['k_{11}(x,' num2str(y1) ')'])
%         %ylim([0.1 0.24])
%     elseif i == 2
%         xlabel('x')
%         ylabel(['k_{11}(x,' num2str(0.5) ')'])
%         ylim([0.11 0.26])
%     elseif i == 3
%         xlabel('x')
%         ylabel(['k_{11}(x,' num2str(0.8) ')'])
%         ylim([0.14 0.28])
%     elseif i == 4
%         xlabel('x')
%         ylabel(['k_{11}(x,' num2str(y2) ')'])
%     end
%     legend('Exact','NL = 2%','NL = 4%','NL = 6%','location','northwest')
%     ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
%     ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);
% end





% K11 PLOTS (VERSION 2)
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
        K = K + Ksol(i,:,j)';
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
    if i == 1 % NL = 0.001
        PlotOptions.LineStyle = '-';
        PlotOptions.Color = 'b';
        PlotOptions.Marker = 'o';
        PlotOptions.MarkerSize = 6;
    elseif i == 2 % NL = 0.02
        PlotOptions.LineStyle = '-';
        PlotOptions.Color = 'k';
        PlotOptions.Marker = '.';
        PlotOptions.MarkerSize = 14;
    elseif i == 3 % NL = 0.05
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
        legend('Exact','NL=0.1%','NL=1.5%','NL=3%','location','northwest')
        xlabel('x')
        title(['k_{11}(x,' num2str(0.2) ')'])
        ylim([0.08 0.26])
    elseif i == 2
        xlabel('x')
        title(['k_{11}(x,' num2str(0.4) ')'])
        ylim([0.10 0.25])
    elseif i == 3
        xlabel('x')
        title(['k_{11}(x,' num2str(0.6) ')'])
        ylim([0.10, 0.25])
        %yticks([0.12, 0.15, 0.2, 0.25])
    elseif i == 4
        xlabel('x')
        title(['k_{11}(x,' num2str(0.8) ')'])
        ylim([0.12 0.25])
        yticks([0.12, 0.15, 0.2, 0.25])
    end
    %legend('Exact','NL=2%','NL=4%','NL=6%','location','northwest')
    ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
    ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);
end



% % K22 PLOTS (not used)
% % plot exact
% PlotOptions.LineWidth = 1.2;
% PlotOptions.LineStyle = '-';
% PlotOptions.Color = 'r';
% PlotOptions.Marker = 'none';
% subplot(2,3,4)
% plot(x,Ke22(0.2*ones(n+1,1),x),PlotOptions), hold on
% subplot(2,3,5)
% plot(x,Ke22(0.5*ones(n+1,1),x),PlotOptions), hold on
% subplot(2,3,6)
% plot(x,Ke22(0.8*ones(n+1,1),x),PlotOptions), hold on
% 
% for i = 1:3
%     K = zeros(2*(n+1)^2,1);
%     for j = 1:repetition
%         K = K + Ksol(i,:,j)';
%     end
%     K = K/repetition; % average
%     K11 = K(1:(n+1)^2);
%     K22 = K((n+1)^2+1:end);
%     K11 = reshape(K11,n+1,n+1);
%     K22 = reshape(K22,n+1,n+1);
%     K11 = griddedInterpolant(X,Y,K11); % to guess values out of cheby's grid
%     K22 = griddedInterpolant(X,Y,K22);
%     
%     % plot options
%     PlotOptions.LineWidth = 1.2;
%     if i == 1 % NL = 0.02
%         PlotOptions.LineStyle = '-';
%         PlotOptions.Color = 'k';
%         PlotOptions.Marker = '.';
%         PlotOptions.MarkerSize = 10;
%     elseif i == 2 % NL = 0.04
%         PlotOptions.LineStyle = '-';
%         PlotOptions.Color = 'b';
%         PlotOptions.Marker = 'o';
%         PlotOptions.MarkerSize = 6;
%     elseif i == 3 % NL = 0.06
%         PlotOptions.LineStyle = '-';
%         PlotOptions.Color = [0 0.8 0.5];
%         PlotOptions.Marker = '*';
%         PlotOptions.MarkerSize = 6;
%     end
% 
%     % plots
%     subplot(2,3,4)
%     %plot(x,K11(x,yline1),PlotOptions), hold on
%     plot(x,K22(0.2*ones(n+1,1),x),PlotOptions), hold on
%     subplot(2,3,5)
%     %plot(x,K22(x,yline1),PlotOptions), hold on
%     plot(x,K22(0.5*ones(n+1,1),x),PlotOptions), hold on
%     %plot(x,K22(x,yline1),PlotOptions), hold on
%     subplot(2,3,6)
%     %plot(x,K11(x,yline2),PlotOptions), hold on
%     plot(x,K22(0.8*ones(n+1,1),x),PlotOptions), hold on
%     
%     
%     %subplot(2,4,4)
%     %plot(x,K22(x,yline2),PlotOptions), hold on
%     %plot(x,K11(x,0.8*ones(n+1,1)),PlotOptions), hold on
%     
% end
% 
% % adjust plots
% for i = 1:3
%     subplot(2,3,i+3)
%     axis square
%     if i == 1
%         legend('Exact','NL = 2%','NL = 4%','NL = 6%','location','northwest')
%         xlabel('x')
%         ylabel(['k_{22}(x,' num2str(y1) ')'])
%         %ylim([0.1 0.24])
%     elseif i == 2
%         xlabel('x')
%         ylabel(['k_{22}(x,' num2str(0.5) ')'])
%         %ylim([0.11 0.26])
%     elseif i == 3
%         xlabel('x')
%         ylabel(['k_{22}(x,' num2str(0.8) ')'])
%         %ylim([0.14 0.28])
%     end
%     legend('Exact','NL = 2%','NL = 4%','NL = 6%','location','northwest')
%     ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
%     ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);
% end



%% CONTOUR PLOTS
figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;

Ke11 = k11(X,Y); % exact K11
Ke22 = k22(X,Y); % exact K22

subplot(2,4,1)
contour(X,Y,Ke11,20)
%surf(X,Y,Ke11), zlim([0, 0.3])
axis square
xlabel('x')
ylabel('y')
yticks([0 0.5 1])
title('Exact')
% subplot(2,4,5)
% contour(X,Y,Ke22,20)
% axis square
% yticks([0 0.5 1])
% title('Exact')

for i = 1:3
    K = zeros(2*(n+1)^2,1);
    for j = 1:repetition
        K = K + Ksol(3*(i-1)+1,:,j)';
    end
    K = K/repetition; % average
    K11 = K(1:(n+1)^2);
    K22 = K((n+1)^2+1:end);
    K11 = reshape(K11,n+1,n+1);
    K22 = reshape(K22,n+1,n+1);
    
    % plots
    subplot(2,4,i+1)
    contour(X,Y,K11,20)
    %surf(X,Y,K11), zlim([0, 0.3])
    axis square
    xlabel('x')
    ylabel('y')
    yticks([0 0.5 1])
    if i == 1
        title('Full grid')
    elseif i == 2
        title('MS1')
    elseif i == 3
        title('MS2')
    end
    
end

ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);



figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;

Ke11 = k11(X,Y); % exact K11
Ke22 = k22(X,Y); % exact K22

subplot(2,4,5)
contour(X,Y,Ke22,20)
axis square
xlabel('x')
ylabel('y')
yticks([0 0.5 1])
title('Exact')

for i = 1:3
    K = zeros(2*(n+1)^2,1);
    for j = 1:repetition
        K = K + Ksol(3*(i-1)+1,:,j)';
    end
    K = K/repetition; % average
    K11 = K(1:(n+1)^2);
    K22 = K((n+1)^2+1:end);
    K11 = reshape(K11,n+1,n+1);
    K22 = reshape(K22,n+1,n+1);
    
    % plots
    subplot(2,4,i+1+4)
    contour(X,Y,K22,20)
    axis square
    xlabel('x')
    ylabel('y')
    yticks([0 0.5 1])
    if i == 1
        title('Full grid')
        hold on
    elseif i == 2
        title('MS1')
%         hold on
%         rectangle('position',[x1 y1 x2-x1 y2-y1],'facecolor',[.8 .8 .8],'edgecolor',[.8 .8 .8])
%         
%         % inner points
%         xi = zeros(n+1,1);
%         yi = zeros(n+1,1);
%         for i = 1:n+1
%             if l1*x(i) >= x1 && l1*x(i) <= x2
%                 xi(i) = 1;
%             end
%             if l2*y(i) >= y1 && l2*y(i) <= y2
%                 yi(i) = 1;
%             end
%         end
% 
%         ell_inner = zeros((n+1)^2,1); % to select inner Jacobian lines
%         for i = 1:n+1
%             for j = 1:n+1
%                 if xi(j) == 1 && yi(i) == 1 
%                     ell = (j-1) + (i-1)*(n+1) + 1;
%                     ell_inner(ell) = ell;
%                 end
%             end
%         end
%         ell_inner = nonzeros(ell_inner);
%         
%         plot(xx(ell_inner),yy(ell_inner),'ko','MarkerSize',3)
    elseif i == 3
        title('MS2')
%         hold on
%         
%         % selection points
%         ell_inner = [51,54,57,60,84,87,90,93,117,120,123,126,147,150,153,156,180,183,186,189]';
%         
%         plot(xx(ell_inner),yy(ell_inner),'k*')
    end
    
end

ht = findall(gcf,'type','text'); set(ht,'fontsize',15);
ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);

































