% run startup
% define parameters for generating synthetic data

par1(1) = 300;
par1(2) = 1;
par1(3) = 0.07;%0.03;
par1(4) = 0.07;%0.03;
par1(5) = 1;
par1(6) = 1;
par1(7) = 0;
 
 
par2(1) = 100;
par2(2) = 3;
par2(3) = 0.03;
par2(4) = 0.03;
par2(5) = 1;
par2(6) = 1;
par2(7) = 18;
 
%%
 
Tfinal = 1500; % 25 hours
Noise1 = sqrt(0.1); 
CellNum = 1000; % number of osc/non-osc cells - total = 2*CellNum
 
[ x, dataNORMED ] = GetTimeSeries( par1,par2, Tfinal, Noise1, CellNum );
 
 
%%
% add SE to data
 
 
% generate SE long term trend...
 
x = (0:30:Tfinal)/60; % in hours
 
% define parameters for generating time series
% note all s.t.d kept at 1
a = exp(-4);%0.05; %0.1
% b = par(2);
cov1 = exp(-a*x.^2);
 
CovMatrix1 = zeros(length(x),length(x));
 
for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);
 
end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';
 
% generate 2*CellNum trends...
MU = zeros(2*CellNum,length(x));
SIGMA = CVM1;
data1 = mvnrnd(MU,SIGMA);

SE = data1';
 
%%
% %add together and check
y = dataNORMED + SE;
time = x';
% plot(x,y(:,cell))
% hold on
% plot(x,dataNORMED(:,cell))
% plot(x,SE(:,cell))
% hold off
 
%%
% for trend of -6
currentdata = y;
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1);
 
parfor i = 1:size(currentdata,2);
    %i
    y1 = currentdata(:,i);  
    x = time;
    x(y1==0) = []; %deletes times from which no signal
    samp = length(x);
    y1(y1==0) = [];
    y1 = y1 - mean(y1);
    
% remove trend from data
 
    Noise = Noise1
%     Noise = Noise1/std(y1);
%     y1 = y1/std(y1);
    raw = y1;
%IMPORTANT PARAMETER - number (3rd input) controls how slow
% trend is for cell data
    [m,par0] = detrenddataNEW(raw,x,-6);
    y1 = y1-m; %detrended y1  
    y1 = y1-mean(y1);
    

% fit OU and OUoscillatory models
 
    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise1);
    par0M(i,:) = par0;
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
 
end
 x = time;
 BICdiffTOT6 = [BICdiffM];
 par0TOT6 = [par0M];
 par1TOT6 = [par1M];
 par2TOT6 = [par2M];
 6
 %%
 % for -4 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
 
parfor i = 1:size(currentdata,2);
%     i
    y1 = currentdata(:,i);  
    x = time;
    x(y1==0) = []; %deletes times from which no signal
    samp = length(x);
    y1(y1==0) = [];
    y1 = y1 - mean(y1);
    
% remove trend from data
    Noise = Noise1
%     Noise = Noise1/std(y1);
%     y1 = y1/std(y1);
    raw = y1;
%IMPORTANT PARAMETER - number (3rd input) controls how slow
% trend is for cell data
    [m,par0] = detrenddataNEW(raw,x,-4);
    y1 = y1-m; %detrended y1  
    y1 = y1-mean(y1);
     
 
% fit OU and OUoscillatory models
 
    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise1);
    par0M(i,:) = par0;
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
 
end
 x = time;
 BICdiffTOT4 = [BICdiffM];
 par0TOT4 = [par0M];
 par1TOT4 = [par1M];
 par2TOT4 = [par2M];
 4
%% 
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
 
parfor i = 1:size(currentdata,2);
%     i
    y1 = currentdata(:,i);  
    x = time;
    x(y1==0) = []; %deletes times from which no signal
    samp = length(x);
    y1(y1==0) = [];
    y1 = y1 - mean(y1);
    
% remove trend from data
 
    Noise = Noise1;
%     Noise = Noise1/std(y1);
%     y1 = y1/std(y1);
    raw = y1;
%IMPORTANT PARAMETER - number (3rd input) controls how slow
% trend is for cell data
    [m,par0] = detrenddataNEW(raw,x,-2);
    y1 = y1-m; %detrended y1  
    y1 = y1-mean(y1);
    
 
% fit OU and OUoscillatory models
 
    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise1);
    par0M(i,:) = par0;
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
 
end
 x = time;
 BICdiffTOT2 = [BICdiffM];
 par0TOT2 = [par0M];
 par1TOT2 = [par1M];
 par2TOT2 = [par2M];
 2
 
%%
 repeats = 1;
 
[ synthOUhier6 ] = MakesynthOUHIERACHICALwithSE(par0TOT6, par1TOT6,repeats,x );
[ BICdiffsynthTOT6 ] = BICdistDATAnewNORENORM( synthOUhier6,time,par1TOT6,repeats, -6);
6
 
[ synthOUhier4 ] = MakesynthOUHIERACHICALwithSE(par0TOT4, par1TOT4,repeats,x );
[ BICdiffsynthTOT4 ] = BICdistDATAnewNORENORM( synthOUhier4,time,par1TOT4,repeats, -4);
4
 
[ synthOUhier2 ] = MakesynthOUHIERACHICALwithSE(par0TOT2, par1TOT2,repeats,x );
[ BICdiffsynthTOT2 ] = BICdistDATAnewNORENORM( synthOUhier2,time,par1TOT2,repeats, -2);
2
%%
disp('done Gillespiewithtrend')
% save('Gillespiewithtrend642rerun.mat') 
% clear


%%
% load data used in paper

% load Fig6


%%
% - caluculate FDRs then plot

% BICdiffTOT = (BICdiffTOT1);
BICdiffTOT = (BICdiffTOT4);
BICdiffsynthTOT = (BICdiffsynthTOT4);

BICdiffTOT(BICdiffTOT<0) = 0;
BICdiffsynthTOT(BICdiffsynthTOT<0) = 0;

upper = max([BICdiffTOT;BICdiffsynthTOT]);
lower1 = min([BICdiffTOT;BICdiffsynthTOT]);
lower = upper - 0.9*(upper-lower1);



range = linspace(lower,upper,20);

for i= 1:length(range)
%     i
    cutoff = range(i);
    num = sum(BICdiffTOT<cutoff)/length(BICdiffTOT);
    denom = sum(BICdiffsynthTOT<cutoff)/length(BICdiffsynthTOT);
    piest(i) =  num/denom;
end
    

figure()
subplot(1,2,1)
plot(range,piest)

% cubic spline regression
xx = linspace(lower,upper,100);
yy = spline(range,piest,xx);
% plot(xx,yy)

% figure()

% hold on
% plot(xx,yy,'color','r')
% ylim([0 1])
% hold off
% xlabel('\lambda')
% ylabel('\pi_0(\lambda)')
% xlim([lower1,upper])
% a = xlim();
% b = ylim();
% text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'a)'},...
%     'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

piGUESS1 = yy(1);


% Go through BIClist calculating q values
% choose initial cutoff level...

[BICdiffTOT,I] = sort(BICdiffTOT);

q1 = zeros(length(BICdiffTOT),1);
% find BIC list location above thresh...

for i = 1:length(BICdiffTOT)

    Thresh = BICdiffTOT(i);
    (sum(BICdiffsynthTOT>Thresh)/length(BICdiffsynthTOT));
    (sum(BICdiffTOT>Thresh)/length(BICdiffTOT));
    q1(i) = piGUESS1*(sum(BICdiffsynthTOT>Thresh)/length(BICdiffsynthTOT))/(sum(BICdiffTOT>Thresh)/length(BICdiffTOT));
end

q = 0.05; % q-value cut-off
cutoff = find(q1<q,1,'first')
[w,l] = sort(I);
Reorderedq = q1(l); % ordered list of q-values
PassList = Reorderedq<q;
FP = sum(PassList(1:1000))  % number of false positives
TP = sum(PassList(1001:2000)) % number of true positives
FDR = FP/(FP+TP) % false discovery rate

 
%%
cell = 1015;
% plot the detrended data
y1 = y(:,cell);
y1 = y1/std(y1);
raw = y1;
%IMPORTANT PARAMETER - number (3rd input) controls how slow
% trend is for cell data
[m1,par0] = detrenddataNEW(raw,x,-6);
det1 = y1-m1; %detrended y1   
 
[m2,par0] = detrenddataNEW(raw,x,-4);
det2 = y1-m2; %detrended y1   
 
[m3,par0] = detrenddataNEW(raw,x,-2);
det3 = y1-m3; %detrended y1   
 
 
%% put together
 
 
subplot(3,3,1)
plot(x,y(:,cell))
hold on
plot(x,dataNORMED(:,cell),'color',[0 0.7 0],'LineStyle','--')
plot(x,SE(:,cell),'color','r')
hold off
xlim([ 0 max(x)])
xlabel('Time (hours)')
ylabel('Normalised expression')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
 
subplot(3,3,2)
plot(x,raw,x,m1)
xlim([ 0 max(x)])
xlabel('Time (hours)')
ylabel('Normalised expression')
t = title('Detrend \alpha_{SE} = exp(-6)')
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(3,3,5)
plot(x,raw,x,m2)
xlim([ 0 max(x)])
xlabel('Time (hours)')
ylabel('Normalised expression')
t = title('Detrend \alpha_{SE} = exp(-4)')
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'D'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(3,3,8)
plot(x,raw,x,m3)
xlim([ 0 max(x)])
xlabel('Time (hours)')
ylabel('Normalised expression')
t = title('Detrend \alpha_{SE} = exp(-2)')
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'F'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
 
 
subplot(3,3,3)
pers1 =2*pi./par2TOT6(CellNum+1:end,2);
av1 = mean(pers1(pers1<6));
edges = linspace(0,6,40);
hist(pers1(pers1<6),edges)
xlabel('Period (hours)')
ylabel('Frequency')
xlim([0 5])
t = title('Estimated period');
t.FontWeight = 'normal';
ylim([0 400])
text(2.5,350,['mean = ',num2str(av1,3)])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(3,3,6)
pers1 =2*pi./par2TOT4(CellNum+1:end,2);
av2 = mean(pers1(pers1<6));
edges = linspace(0,6,40);
hist(pers1(pers1<6),edges)
xlabel('Period (hours)')
ylabel('Frequency')
xlim([0 5])
ylim([0 400])
text(2.5,350,['mean = ',num2str(av2,3)])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'E'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(3,3,9)
pers1 =2*pi./par2TOT2(CellNum+1:end,2);
av3 = mean(pers1(pers1<6));
edges = linspace(0,6,40);
hist(pers1(pers1<6),edges)
xlabel('Period (hours)')
ylabel('Frequency')
xlim([0 5])
ylim([0 400])
text(2.5,350,['mean = ',num2str(av3,3)])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'G'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')