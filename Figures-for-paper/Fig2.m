% generate an OUosc time series...

x = linspace(0,100,100);
% define parameters for generating time series
% note all s.t.d kept at 1
par(1) = 0.15;
par(2) = 0.5;

a = par(1);%0.05; %0.1
b = par(2);
Q2 = b/(2*pi()*a);
% cov1 = exp(-a*x);
cov1 = exp(-a*x).*cos(b*x);
% cov2 = exp(-0.05*x.^2);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


MU = zeros(1,length(x));
% Noise parameter represents standard dev of meas noise
Noise = 0.1;
Meas = diag((Noise^2).*ones(1,length(x)));
CVM1 = CVM1 + Meas;
SIGMA = CVM1; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);
OUosc = data1;

%%
% get BICdiff
[BICdiff1, par1a, par2a] = getBICdiff(x',OUosc',Noise);

% x = x';
% y1 = y1';

plot(x,data1)

%%

% generate SE long term trend...

x = linspace(0,100,100);
% define parameters for generating time series
% note all s.t.d kept at 1
par(1) = 0.001;
par(2) = 0.5;

a = par(1);%0.05; %0.1
% b = par(2);
% Q2 = b/(2*pi()*a);
% cov1 = exp(-a*x);
% cov1 = exp(-a*x).*cos(b*x);
cov1 = exp(-a*x.^2);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';

% generate 3 examples...
MU = zeros(3,length(x));
SIGMA = CVM1; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);


% x = x';
% y1 = y1';
SE = data1;
plot(x,SE(1,:))
hold on
plot(x,SE(2,:),':','color','r')
plot(x,SE(3,:),':','color','r')
hold off

%%
% plot two time series added together, and calculate BIC

SigAndTrend = OUosc + SE(1,:);

[BICdiff2, par1b, par2b] = getBICdiff(x',SigAndTrend',Noise);


%%
plot(x,SigAndTrend)


%%
% detrend data and recalculate BICdiff

[m] = detrenddata(SigAndTrend',x',-5);
y1 = SigAndTrend-m'; %detrended y1  
plot(x,y1)

[BICdiff3, par1c, par2c] = getBICdiff(x',y1',Noise);

%%
% load fig2
% plot everything

subplot(2,2,1)
plot(x,OUosc)
xlabel('Time (hours)')
ylabel('Original signal')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
str1 = sprintf('%.2f',par2a(2));
str2 = sprintf('%.2f',BICdiff1);
t = title(['\beta = ',str1,', LLR = ',str2]);
t.FontWeight = 'normal';

subplot(2,2,2)
plot(x,SE(1,:),'color','k')
hold on
plot(x,SE(2,:),'--','color','r')
plot(x,SE(3,:),'--','color','r')
hold off
xlabel('Time (hours)')
ylabel('Long term trend')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


subplot(2,2,3)
plot(x,SigAndTrend)
hold on
plot(x,SE(1,:),'--','color','k')
plot(x,m,'--','color',[0 0.8 0 ])
hold off
xlabel('Time (hours)')
ylabel('Signal and trend')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
str1 = sprintf('%.2f',par2b(2));
str2 = sprintf('%.2f',BICdiff2);
t = title(['\beta = ',str1,', LLR = ',str2]);
t.FontWeight = 'normal';

subplot(2,2,4)
plot(x,y1)
xlabel('Time (hours)')
ylabel('Detrended signal')
ylim([ -3 4])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'D'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
str1 = sprintf('%.2f',par2c(2));
str2 = sprintf('%.2f',BICdiff3);
t = title(['\beta = ',str1,', LLR = ',str2]);
t.FontWeight = 'normal';