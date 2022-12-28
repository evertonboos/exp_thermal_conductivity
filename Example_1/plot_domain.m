% Routine to plot domain and selected points for figure (example) in the
% beginning of section 3.2 (or close), to enforce possible choices of
% measurement points given to the code.

clear; close all;

% plot positioning
aa = (1920-1366)/1920*0.5; % notebook -> monitor display
bb = (1080-768)/1080*0.5;
%posi = [0 0 1 1]; % full screen (any monitor)
posi = [aa bb 1366/1920 768/1080]; % centered, notebook size

%% Parameters
n = 9;
[~,x] = cheby1(n,0,1); % Cheby's points
l1 = 1; % x limit
l2 = 1; % y limit

% grid
[X,Y] = ndgrid(x,x);
xx = X(:);
yy = Y(:);

%% Plot Inner domain
figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;

% inner region limits
x1 = 0.08;
x2 = 0.92;
y1 = 0.2;
y2 = 0.8;

% plots
subplot(1,3,2)
rectangle('position',[0 0 l1 l2],'linewidth',1.5) % general domain
hold on
rectangle('position',[x1 y1 x2-x1 y2-y1],'facecolor',[.8 .8 .8],'edgecolor',[.8 .8 .8])
plot(X,Y,'kx')

indexes = [33,34,35,36,37,38,43,44,45,46,47,48,...
    53,54,55,56,57,58,63,64,65,66,67,68];
points = {'33','34','35','36','37','38',...
    '43','44','45','46','47','48',...
    '53','54','55','56','57','58',...
    '63','64','65','66','67','68'};
plot(xx(indexes),yy(indexes),'k.','markersize',20)
for i = 1:length(indexes)
    text(xx(indexes(i))-0.04, yy(indexes(i))+0.06, points{i})
end


xlim([-0.1 1.1])
ylim([-0.1 1.1])
axis square
xticks([0 0.5 1])
yticks([0 0.5 1])
xlabel('x axis')
ylabel('y axis')

ht = findall(gcf,'type','text'); set(ht,'fontsize',12);
ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);

%% Plot Selection
figure,
fig = gcf;
fig.Units = 'normalized';
fig.OuterPosition = posi;

% plots
subplot(1,3,2)
rectangle('position',[0 0 l1 l2],'linewidth',1.5) % general domain
hold on
plot(X,Y,'kx')

indexes = [33,35,37,44,46,48,53,55,57,64,66,68];
points = {'33','35','37','44','46','48','53','55','57','64','66','68'};
plot(xx(indexes),yy(indexes),'k.','markersize',20)
for i = 1:length(indexes)
    text(xx(indexes(i))-0.04, yy(indexes(i))+0.06, points{i})
end


xlim([-0.1 1.1])
ylim([-0.1 1.1])
axis square
xticks([0 0.5 1])
yticks([0 0.5 1])
xlabel('x axis')
ylabel('y axis')

ht = findall(gcf,'type','text'); set(ht,'fontsize',12);
ht = findall(gcf,'type','axes'); set(ht,'fontsize',14);





















