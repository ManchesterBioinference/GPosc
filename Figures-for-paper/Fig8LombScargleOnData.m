%loads data from excel spreadsheet - need to change file name
num = xlsread('Hes1ublucF5_sparse_combine_30minexp_forBayesian.xlsx',1);
num(isnan(num)) = 0;
 
 
% assigns time 
time = num(:,1);
time = time/(60*60*1000); %convert from ms to hours
time1 = time; % store for later...
 
 
 
% Manually enter the columns of data and bckgd for all experiments
data.exp1.data = num(:,2:8); data.exp1.bckgd = num(:,9:12);
 
data.exp2.data = num(:,13:24); data.exp2.bckgd = num(:,25:28);
 
 
% creates empty vectors where final results will be stored
BICdiffTOT1 = [];
 
 
%%
% gets number of experiments
fields = fieldnames(data)
 
%%
% creat BICdiffTOT1, par1TOT1 and par2TOT1 for first (C17)
% data set
tic
for j =1:length(fields)
j
current = data.(fields{j})
 
 
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
 
    y1 = y1/std(y1);
    
%IMPORTANT PARAMETER - number (3rd input) controls how slow
% trend is for cell data
    raw = y1;
    [m] = detrenddata(raw,x,-4.5);
    y1 = y1-m; %detrended y1   
%     y1 = y1/std(y1);
    
    [f,P,prob] = lomb(x,y1,3,1.05); 
 
    BICdiffM(i) = min(prob);
    
% show figure - comment out "showfigure" is don't want 
%individual cell figures popping up
%     showfigure(x,m,raw,y1,BICdiff,par1,par2,i)
 
end
% updates total list for all pooled experiments
 BICdiffTOT1 = [BICdiffTOT1;BICdiffM];
 
end
 
toc

%%
 
num = xlsread('Controlpromoter_cellsforBayesian_030615.xlsx',1);
num(isnan(num)) = 0;
 
% assigns time 
time = num(:,1);
time = time/(60*60*1000); %convert from ms to hours
%%
 
data.exp1.data = num(:,2:25); data.exp1.bckgd = num(:,26:29);
 
data.exp2.data = num(:,30); data.exp2.bckgd = num(:,31:34);
 
% data.exp3.data = num(:,29:55); data.exp3.bckgd = num(:,56:59);
 
% data.exp4.data = num(:,23); data.exp4.bckgd = num(:,25:28);
% 
% data.exp5.data = num(:,23); data.exp5.bckgd = num(:,25:28);
 
 
% creates empty vectors where final results will be stored
BICdiffTOT2 = [];
 
 
 
 
% gets number of experiments
fields = fieldnames(data)
for j =1:length(fields)
j
current = data.(fields{j})
 
 
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
 
    y1 = y1/std(y1);
    
%IMPORTANT PARAMETER - number (3rd input) controls how slow
% trend is for cell data
%     raw = y1;
%     [m] = detrenddata(raw,x,-4.5);
%     y1 = y1-m; %detrended y1   
%     y1 = y1/std(y1);
    
       [f,P,prob] = lomb(x,y1,3,1.05); 
 
    BICdiffM(i) = min(prob);
% show figure - comment out "showfigure" is don't want 
%individual cell figures popping up
%     showfigure(x,m,raw,y1,BICdiff,par1,par2,i)
 
end
% updates total list for all pooled experiments
 BICdiffTOT2 = [BICdiffTOT2;BICdiffM];
 
 
end
 
%%
BICdiffTOT = (BICdiffTOT2); % for control data

%%

% With B-Hochberg
% first need ordered p-values

[pvalsO,I] = sort(BICdiffTOT);

q = 0.05; % choose cutoff q-value to control FDR

for k = 1:length(pvalsO)
    pass(k) = pvalsO(k)<(q*k/length(pvalsO));
end

find(pass>0,1,'last')

[w,l] = sort(I);

OrderedPassList = pass(l); % this says for each cell whether it passes


