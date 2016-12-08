x = linspace(0,100,1000);
%%
% define parameters for generating time series
% note all s.t.d kept at 1
par(1) = 0.5;
par(2) = 0.01;
par(3) = 0.04;
par(4) = 0.5;
par(5) = 0.2;
par(6) = 0.5;

%%

a = par(1);%0.05; %0.1

cov1 = exp(-a*x);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


MU = zeros(1,length(x));
SIGMA = CVM1; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);


OU1 = data1;
plot(x,data1)

%%
a = par(2);%0.05; %0.1

cov1 = exp(-a*x);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


MU = zeros(1,length(x));
SIGMA = CVM1; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);

OU2 = data1;
plot(x,data1)

%%

a = par(3);%0.05; %0.1
b = par(4);
Q1 = b/(2*pi()*a);

cov1 = exp(-a*x).*cos(b*x);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


MU = zeros(1,length(x));
SIGMA = CVM1; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);

OUosc1 = data1;
plot(x,data1)

%%
a = par(5);%0.05; %0.1
b = par(6);
Q2 = b/(2*pi()*a);

cov1 = exp(-a*x).*cos(b*x);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


MU = zeros(1,length(x));
SIGMA = CVM1; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);

OUosc2 = data1;
plot(x,data1)
%%
% to reproduce in paper
% load fig1
% 
a = par(3);%0.05; %0.1
b = par(4);
Q1 = b/(2*pi()*a);

a = par(5);%0.05; %0.1
b = par(6);
Q2 = b/(2*pi()*a);

%%

% change sine wave starting poisition to approximate 
deter1 = sqrt(2)*sin(par(4).*x+3.5);
deter2 = sqrt(2)*sin(par(6).*x);
% plot(x,deter)

subplot(2,2,1)
plot(x,OU1)
xlabel('Time (hours)')
ylabel('y(time)')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
str = sprintf('%.2f',par(1));
t = title(['\alpha = ',str]);
t.FontWeight = 'normal';
subplot(2,2,2)
plot(x,OU2)
xlabel('Time (hours)')
ylabel('y(time)')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
str = sprintf('%.2f',par(2));
t = title(['\alpha = ',str]);
t.FontWeight = 'normal';
subplot(2,2,3)
plot(x,OUosc1)
hold on
plot(x,deter1,':','color','r')
hold off
xlabel('Time (hours)')
ylabel('y(time)')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
str1 = sprintf('%.2f',par(3));
str2 = sprintf('%.2f',par(4));
str3 = sprintf('%.2f',Q1);
t = title(['\alpha = ',str1,', \beta = ',str2,', Quality = ',str3]);
t.FontWeight = 'normal';
subplot(2,2,4)
plot(x,OUosc2)
hold on
plot(x,deter2,':','color','r')
hold off
xlabel('Time (hours)')
ylabel('y(time)')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'D'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
str1 = sprintf('%.2f',par(5));
str2 = sprintf('%.2f',par(6));
str3 = sprintf('%.2f',Q2);
t = title(['\alpha = ',str1,', \beta = ',str2,', Quality = ',str3]);
t.FontWeight = 'normal';


