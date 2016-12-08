% set (parcluster, 'NumWorkers', 12)
% matlabpool open 12
%%
% generate synthetic data using OU GP model only
CellNum = 1000;
Tfinal = 25*60;
time = (0:30:Tfinal)/60;
Noise1 = sqrt(0.1);

for j = 1:CellNum

x = time;
v = 1;%par1TOT(2);
a = 1;%par1TOT(1);
cov1 = v*exp(-a*x);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';

% generate 3 examples...
MU = zeros(1,length(x));
Meas = diag((Noise1^2).*ones(1,length(x)));
SIGMA = CVM1+Meas; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);


% x = x';
% y1 = y1';
Synthcurr = data1';
Synth1(:,j) =Synthcurr;

end


Synth = [Synth1];%,Synth2];
% x = time';
%%
% run through pipeline

y = Synth;


currentdata = y;
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1);


% load cell data for current experiment - loops through cells
parfor i = 1:size(currentdata,2);
%     i
    y1 = currentdata(:,i);  
     x = time';
    x(y1==0) = []; %deletes times from which no signal
    samp = length(x);
    y1(y1==0) = [];
    y1 = y1 - mean(y1);

    Noise = Noise1/std(y1);
    y1 = y1/std(y1);

% fit OU and OUoscillatory models

    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise);
%     par0M(i,:) = par0;
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;


end

%%
 x = time';
 BICdiffTOT = [BICdiffM];
%  par0TOT = [par0M];
 par1TOT = [par1M];
 par2TOT = [par2M];
 
%  tic
 
repeats = 1;
% turn on parallel
% [ synthOUhier ] = MakesynthOUHIERACHICALwithSE(par0TOT, par1TOT,repeats,x );
[ synthOUhier ] = MakesynthOUHIERACHICAL(par1TOT,repeats,x );
% [ BICdiffsynthTOT ] = BICdistDATAnew( synthOUhier,time,par1TOT,repeats, -7);
[ BICdiffsynthTOT ] = BICdistDATA( synthOUhier,x,par1TOT,repeats);
% toc
disp('done OUnotrendshort')
% save('OUnotrendshort.mat') 
clear
 %%
% do again for longer time (50 hours)
 CellNum = 1000;
Tfinal = 50*60; %50 hours
time = (0:30:Tfinal)/60;
Noise1 = sqrt(0.1);

for j = 1:CellNum

x = time;
v = 1;%par1TOT(2);
a = 1;%par1TOT(1);
cov1 = v*exp(-a*x);
CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';

% generate 3 examples...
MU = zeros(1,length(x));
Meas = diag((Noise1^2).*ones(1,length(x)));
SIGMA = CVM1+Meas; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);


% x = x';
% y1 = y1';
Synthcurr = data1';
Synth1(:,j) =Synthcurr;

end


Synth = [Synth1];%,Synth2];
% x = time';
%%
% run through pipeline

y = Synth;


currentdata = y;
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1);


% load cell data for current experiment - loops through cells
parfor i = 1:size(currentdata,2);
%     i
    y1 = currentdata(:,i);  
     x = time';
    x(y1==0) = []; %deletes times from which no signal
    samp = length(x);
    y1(y1==0) = [];
    y1 = y1 - mean(y1);
    
    Noise = Noise1/std(y1);
    y1 = y1/std(y1);
%     raw = y1;

% fit OU and OUoscillatory models

    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise);
%     par0M(i,:) = par0;
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;


end

%%
 x = time';
 BICdiffTOT = [BICdiffM];
%  par0TOT = [par0M];
 par1TOT = [par1M];
 par2TOT = [par2M];
 
%  tic
 
repeats = 1;
% turn on parallel
% [ synthOUhier ] = MakesynthOUHIERACHICALwithSE(par0TOT, par1TOT,repeats,x );
[ synthOUhier ] = MakesynthOUHIERACHICAL(par1TOT,repeats,x );
% [ BICdiffsynthTOT ] = BICdistDATAnew( synthOUhier,time,par1TOT,repeats, -7);
[ BICdiffsynthTOT ] = BICdistDATA( synthOUhier,x,par1TOT,repeats);
% toc
disp('done OUnotrendlong')
% save('OUnotrendlong.mat') 
clear

%%

% for OU with no trend

load FigS4s % short data

BICdiffTOT(BICdiffTOT<0) = 0;
BICdiffsynthTOT(BICdiffsynthTOT<0) = 0;

subplot(4,2,1)
hist(BICdiffTOT(1:1000),0:30)
t = title(['mean = ',num2str(mean(BICdiffTOT(1:1000)),'%.2f')]);
t.FontWeight = 'normal';
xlim([0 30])
ylim([0 700])
xlabel('LLR')
ylabel('Frequency')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(4,2,3)
hist(BICdiffsynthTOT(1:1000),0:30)
xlim([0 30])
ylim([0 700])
xlabel('LLR')
ylabel('Frequency')
t = title(['mean = ',num2str(mean(BICdiffsynthTOT(1:1000)),'%.2f')]);
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(4,2,5)
qqplot((BICdiffsynthTOT(1:1000)),(BICdiffTOT(1:1000)))
xlabel('Quantiles bootstrap OU data LLRs')
ylabel('Quantiles Gillespie LLRs')
xlim([0 34])
ylim([0 34])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'E'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


subplot(4,2,7)
hist(par1TOT(:,1),0:0.5:5)
t = title(['\alpha estimate (true=1). Mean = ',num2str(mean(par1TOT(:,1)),'%.2f')]);
t.FontWeight = 'normal';
xlabel('\alpha')
ylabel('Frequency')
xlim([0 5])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'G'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


load FigS4l % long data

subplot(4,2,2)
hist(BICdiffTOT(1:1000),0:30)
t = title(['mean = ',num2str(mean(BICdiffTOT(1:1000)),'%.2f')]);
t.FontWeight = 'normal';
xlim([0 30])
ylim([0 500])
xlabel('LLR')
ylabel('Frequency')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(4,2,4)
hist(BICdiffsynthTOT(1:1000),0:30)
xlim([0 30])
ylim([0 500])
xlabel('LLR')
ylabel('Frequency')
t = title(['mean = ',num2str(mean(BICdiffsynthTOT(1:1000)),'%.2f')]);
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'D'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(4,2,6)
qqplot((BICdiffsynthTOT(1:1000)),(BICdiffTOT(1:1000)))
xlabel('Quantiles bootstrap OU data LLRs')
ylabel('Quantiles Gillespie LLRs')
xlim([0 23])
ylim([0 23])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'F'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


subplot(4,2,8)
hist(par1TOT(:,1),0:0.5:5)
t = title(['\alpha estimate (true=1). Mean = ',num2str(mean(par1TOT(:,1)),'%.2f')]);
t.FontWeight = 'normal';
xlabel('\alpha')
ylabel('Frequency')
xlim([0 5])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.06*(b(2)-b(1)),{'H'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')