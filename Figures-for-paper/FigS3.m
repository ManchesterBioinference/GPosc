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
CellNum = 2000; % number of osc/non-osc cells - total = 2*CellNum
 
[ x, dataNORMED ] = GetTimeSeries( par1,par2, Tfinal, Noise1, CellNum );
 
 
%%
% %add together and check
y = dataNORMED;% + SE;
% time = x';
time = x;
% plot(x,y(:,cell))
% hold on
% plot(x,dataNORMED(:,cell))
% plot(x,SE(:,cell))
% hold off

%%
% create a list of potential detrending parameters to scan over
Period = 2.5;
Trendlist = [0.3 0.5 1 2 3 4 5 6]*Period;
aList = log(1./(2*Trendlist.^2));
 
%%
% for trend of -6
currentdata = y;
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1);
currdetrend = aList(1);
% tic  
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
    y1 = y1-m; %detrended y1  
    y1 = y1-mean(y1);
    

% fit OU and OUoscillatory models
 
    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise1);
    par0M(i,:) = par0;
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
 
end
% toc
 x = time;
 BICdiffTOT1 = [BICdiffM];
 par0TOT1 = [par0M];
 par1TOT1 = [par1M];
 par2TOT1 = [par2M];
 1
 %%
 % for -4 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(2);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(3);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 BICdiffTOT3 = [BICdiffM];
 par0TOT3 = [par0M];
 par1TOT3 = [par1M];
 par2TOT3 = [par2M];
 3
 
%%
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(4);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
currdetrend = aList(5);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 BICdiffTOT5 = [BICdiffM];
 par0TOT5 = [par0M];
 par1TOT5 = [par1M];
 par2TOT5 = [par2M];
 5
%%
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(6);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(7);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 BICdiffTOT7 = [BICdiffM];
 par0TOT7 = [par0M];
 par1TOT7 = [par1M];
 par2TOT7 = [par2M];
 7
%% 
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(8);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 BICdiffTOT8 = [BICdiffM];
 par0TOT8 = [par0M];
 par1TOT8 = [par1M];
 par2TOT8 = [par2M];
 8 
%%
 repeats = 1;
 
[ synthOUhier1 ] = MakesynthOUHIERACHICALwithSE(par0TOT1, par1TOT1,repeats,x );
[ BICdiffsynthTOT1 ] = BICdistDATAnewNORENORM( synthOUhier1,time,par1TOT1,repeats, aList(1));
1 
 
[ synthOUhier2 ] = MakesynthOUHIERACHICALwithSE(par0TOT2, par1TOT2,repeats,x );
[ BICdiffsynthTOT2 ] = BICdistDATAnewNORENORM( synthOUhier2,time,par1TOT2,repeats, aList(2));
2 

[ synthOUhier3 ] = MakesynthOUHIERACHICALwithSE(par0TOT3, par1TOT3,repeats,x );
[ BICdiffsynthTOT3 ] = BICdistDATAnewNORENORM( synthOUhier3,time,par1TOT3,repeats, aList(3));
3
 
[ synthOUhier4 ] = MakesynthOUHIERACHICALwithSE(par0TOT4, par1TOT4,repeats,x );
[ BICdiffsynthTOT4 ] = BICdistDATAnewNORENORM( synthOUhier4,time,par1TOT4,repeats, aList(4));
4

[ synthOUhier5 ] = MakesynthOUHIERACHICALwithSE(par0TOT5, par1TOT5,repeats,x );
[ BICdiffsynthTOT5 ] = BICdistDATAnewNORENORM( synthOUhier5,time,par1TOT5,repeats, aList(5));
5

[ synthOUhier6 ] = MakesynthOUHIERACHICALwithSE(par0TOT6, par1TOT6,repeats,x );
[ BICdiffsynthTOT6 ] = BICdistDATAnewNORENORM( synthOUhier6,time,par1TOT6,repeats, aList(6));
6

[ synthOUhier7 ] = MakesynthOUHIERACHICALwithSE(par0TOT7, par1TOT7,repeats,x );
[ BICdiffsynthTOT7 ] = BICdistDATAnewNORENORM( synthOUhier7,time,par1TOT7,repeats, aList(7));
7

[ synthOUhier8 ] = MakesynthOUHIERACHICALwithSE(par0TOT8, par1TOT8,repeats,x );
[ BICdiffsynthTOT8 ] = BICdistDATAnewNORENORM( synthOUhier8,time,par1TOT8,repeats, aList(8));
8
%% 
save('FigS3notrend')

%%
% add SE to data
 
 
% generate SE long term trend...
 
x = (0:30:Tfinal)/60; % in hours
 
% define parameters for generating time series
% note all s.t.d kept at 1
a = exp(-4);%0.05; %0.1
% b = par(2);
cov1 = 2*exp(-a*x.^2);
 
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
currdetrend = aList(1);
tic  
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
    y1 = y1-m; %detrended y1  
    y1 = y1-mean(y1);
    

% fit OU and OUoscillatory models
 
    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise1);
    par0M(i,:) = par0;
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
 
end
toc
 x = time;
 BICdiffTOT1 = [BICdiffM];
 par0TOT1 = [par0M];
 par1TOT1 = [par1M];
 par2TOT1 = [par2M];
 1
 %%
 % for -4 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(2);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(3);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 BICdiffTOT3 = [BICdiffM];
 par0TOT3 = [par0M];
 par1TOT3 = [par1M];
 par2TOT3 = [par2M];
 3
 
%%
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(4);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
currdetrend = aList(5);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 BICdiffTOT5 = [BICdiffM];
 par0TOT5 = [par0M];
 par1TOT5 = [par1M];
 par2TOT5 = [par2M];
 5
%%
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(6);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(7);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 BICdiffTOT7 = [BICdiffM];
 par0TOT7 = [par0M];
 par1TOT7 = [par1M];
 par2TOT7 = [par2M];
 7
%% 
 % for -2 detrend
par0M = zeros(size(currentdata,2),3);
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1); 
currdetrend = aList(8);
 
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
    [m,par0] = detrenddataNEW(raw,x,currdetrend);
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
 BICdiffTOT8 = [BICdiffM];
 par0TOT8 = [par0M];
 par1TOT8 = [par1M];
 par2TOT8 = [par2M];
 8 
%%
 repeats = 1;
 
[ synthOUhier1 ] = MakesynthOUHIERACHICALwithSE(par0TOT1, par1TOT1,repeats,x );
[ BICdiffsynthTOT1 ] = BICdistDATAnewNORENORM( synthOUhier1,time,par1TOT1,repeats, aList(1));
1 
 
[ synthOUhier2 ] = MakesynthOUHIERACHICALwithSE(par0TOT2, par1TOT2,repeats,x );
[ BICdiffsynthTOT2 ] = BICdistDATAnewNORENORM( synthOUhier2,time,par1TOT2,repeats, aList(2));
2 

[ synthOUhier3 ] = MakesynthOUHIERACHICALwithSE(par0TOT3, par1TOT3,repeats,x );
[ BICdiffsynthTOT3 ] = BICdistDATAnewNORENORM( synthOUhier3,time,par1TOT3,repeats, aList(3));
3
 
[ synthOUhier4 ] = MakesynthOUHIERACHICALwithSE(par0TOT4, par1TOT4,repeats,x );
[ BICdiffsynthTOT4 ] = BICdistDATAnewNORENORM( synthOUhier4,time,par1TOT4,repeats, aList(4));
4

[ synthOUhier5 ] = MakesynthOUHIERACHICALwithSE(par0TOT5, par1TOT5,repeats,x );
[ BICdiffsynthTOT5 ] = BICdistDATAnewNORENORM( synthOUhier5,time,par1TOT5,repeats, aList(5));
5

[ synthOUhier6 ] = MakesynthOUHIERACHICALwithSE(par0TOT6, par1TOT6,repeats,x );
[ BICdiffsynthTOT6 ] = BICdistDATAnewNORENORM( synthOUhier6,time,par1TOT6,repeats, aList(6));
6

[ synthOUhier7 ] = MakesynthOUHIERACHICALwithSE(par0TOT7, par1TOT7,repeats,x );
[ BICdiffsynthTOT7 ] = BICdistDATAnewNORENORM( synthOUhier7,time,par1TOT7,repeats, aList(7));
7

[ synthOUhier8 ] = MakesynthOUHIERACHICALwithSE(par0TOT8, par1TOT8,repeats,x );
[ BICdiffsynthTOT8 ] = BICdistDATAnewNORENORM( synthOUhier8,time,par1TOT8,repeats, aList(8));
8
 
%% 
save('FigS3withtrend')
%%
disp('done Gillespiewithtrend')


%%


load('FigS3notrend.mat')
BICdiffTOTm = [BICdiffTOT1,BICdiffTOT2,BICdiffTOT3,BICdiffTOT4,BICdiffTOT5,BICdiffTOT6,BICdiffTOT7,BICdiffTOT8];
BICdiffsynthm = [BICdiffsynthTOT1,BICdiffsynthTOT2,BICdiffsynthTOT3,BICdiffsynthTOT4,BICdiffsynthTOT5,BICdiffsynthTOT6,BICdiffsynthTOT7,BICdiffsynthTOT8];

for i=1:length(aList)
    [ FP(i), TP(i), FDR(i) ] = GetFDRlong( BICdiffTOTm(:,i), BICdiffsynthm(:,i));
end

figure()
subplot(3,2,1)
plot(Trendlist,FP,'x')
title('With no added trend','fontweight','normal')
ylabel('False positive rate')
xlabel('Detrending length scale (hours)')
ylim([0 0.1])
subplot(3,2,3)
plot(Trendlist,TP,'x')
% ylim([900 2000])
ylabel('Statistical power')
xlabel('Detrending length scale (hours)')
subplot(3,2,5)
plot(Trendlist,FDR,'x')
ylim([0 0.1])
ylabel('FDR')
xlabel('Detrending length scale (hours)')
% load other data set
load('FigS3withtrend.mat')
BICdiffTOTm = [BICdiffTOT1,BICdiffTOT2,BICdiffTOT3,BICdiffTOT4,BICdiffTOT5,BICdiffTOT6,BICdiffTOT7,BICdiffTOT8];
BICdiffsynthm = [BICdiffsynthTOT1,BICdiffsynthTOT2,BICdiffsynthTOT3,BICdiffsynthTOT4,BICdiffsynthTOT5,BICdiffsynthTOT6,BICdiffsynthTOT7,BICdiffsynthTOT8];

for i=1:length(aList)
    [ FP(i), TP(i), FDR(i) ] = GetFDRlong( BICdiffTOTm(:,i), BICdiffsynthm(:,i));
end

subplot(3,2,2)
plot(Trendlist,FP,'x')
title('With added trend','fontweight','normal')
ylabel('False positive rate')
xlabel('Detrending length scale (hours)')
ylim([0 0.1])
subplot(3,2,4)
plot(Trendlist,TP,'x')
% ylim([900 2000])
ylabel('Statistical power')
xlabel('Detrending length scale (hours)')
subplot(3,2,6)
plot(Trendlist,FDR,'x')
ylim([0 0.1])
ylabel('FDR')
xlabel('Detrending length scale (hours)')

