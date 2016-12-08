
%%
%loads data from excel spreadsheet - need to change file name to file that
%you want to analyse
num = xlsread('Controlpromoter_longtimecells_forNick.xlsx',1);
num(isnan(num)) = 0;

%%
% assigns time 
time = num(:,1); % enter column of time
time = time/(60*60*1000); %convert from ms to hours
%%
% Manually enter the columns of data and bckgd for all experiments
data.exp1.data = num(:,2:14); data.exp1.bckgd = num(:,17:20);

data.exp2.data = num(:,23); data.exp2.bckgd = num(:,25:28);

% data.exp3.data = num(:,23); data.exp3.bckgd = num(:,25:28);

% data.exp4.data = num(:,23); data.exp4.bckgd = num(:,25:28);
% 
% data.exp5.data = num(:,23); data.exp5.bckgd = num(:,25:28);

%%
% creates empty vectors where final results will be stored
BICdiffTOT = [];
par1TOT = [];
par2TOT = [];

%%
% gets number of experiments
fields = fieldnames(data)

%%

for j =1:length(fields)
j
current = data.(fields{j})

% find background standard deviation
% loads bckgd for current experiment
bckgd = current.bckgd;
% mean of the standard deviation of all 4 cells
% mean(std(bckgd))

% remove slow trend from bckgd 
% number (3rd entry) changes how slow trend line 
% is for BACKGROUND data

[stdev] = bckgdstdev(bckgd,time,-7);

par1M = zeros(size(current.data,2),3);
par2M = zeros(size(current.data,2),4);
BICdiffM = zeros(size(current.data,2),1);

% load cell data for current experiment - loops through cells
for i = 1:size(current.data,2);
    i
    y1 = current.data(:,i);  
    x = time;
    x(y1==0) = []; %deletes times from which no signal
    samp = length(x);
    y1(y1==0) = [];
    y1 = y1 - mean(y1);
    
% remove trend from data

    Noise = stdev/std(y1);
    y1 = y1/std(y1);
    raw = y1;
%IMPORTANT PARAMETER - number (3rd input) controls how slow
% trend is for cell data
    [m] = detrenddata(raw,x,-5);
    y1 = y1-m; %detrended y1   
    y1 = y1/std(y1);        

% fit OU and OUoscillatory models

    [BICdiff, par1, par2] = getBICdiff(x,y1,Noise);
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
% show figure - comment out "showfigure" is don't want 
%individual cell figures popping up
    showfigure(x,m,raw,y1,BICdiff,par1,par2,i)

end
% updates total list for all pooled experiments
 BICdiffTOT = [BICdiffTOT;BICdiffM];
 par1TOT = [par1TOT;par1M];
 par2TOT = [par2TOT;par2M];

end

pass = sum(BICdiffTOT>3)

%%
% plots distribution of BIC scores
figure()
hist(BICdiffTOT)
xlabel('BIC score')
ylabel('frequency')
title('Distribution of BIC scores')
xlim([-5,25])

%%
%only plot periods of those passing
passlist = BICdiffTOT>3;
periods = 2*pi()./par2TOT(:,2);
figure()
hist(periods(passlist),[0:10])
title('Periods of passing cells')
xlabel('period (hours)')
ylabel('frequency')

%%
% plots distribution of quality parameter of all cells
figure()
title('Quality')
hist(1./par2TOT(:,1),[0:0.2:2])
xlim([0 3])
title('Quality of fitted oscillations to all cells')
xlabel('length scale (hours)')
ylabel('frequency')