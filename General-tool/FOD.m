startup
clear
%%

%loads data from excel spreadsheet - need to change file name to file that
%you want to analyse
num = xlsread('Controlpromoter_longtimecells_forNick.xlsx',1);
num(isnan(num)) = 0;


% assigns time 
time = num(:,1); % enter column of time
time = time/(60*60*1000); %convert from ms to hours
time1 = time; % store for later...

% Manually enter the columns of data and associated bckgd for each experiment
data.exp1.data = num(:,2:8); data.exp1.bckgd = num(:,9:12);

data.exp2.data = num(:,13:24); data.exp2.bckgd = num(:,25:28);

% data.exp3.data = num(:,29:55); data.exp3.bckgd = num(:,56:59);

% data.exp4.data = num(:,23); data.exp4.bckgd = num(:,25:28);
% 
% data.exp5.data = num(:,23); data.exp5.bckgd = num(:,25:28);



% creates empty vectors where final results will be stored
BICdiffTOT1 = [];
par0TOT1 = [];
par1TOT1 = [];
par2TOT1 = [];


% gets number of experiments
fields = fieldnames(data)

%%

% choose lengthscale to detrend (in hours)
Lengthscale = 7.5; % roughly 3x period expected (in hours)
DetrendParam = log(1./(2*Lengthscale.^2));

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

par0M = zeros(size(current.data,2),3); %sets up matrices to fill
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
%     Noise = stdev;
    Noise = stdev/std(y1);
    y1 = y1/std(y1);
    raw = y1;
%IMPORTANT PARAMETER - number (3rd input) controls how slow
% trend is for cell data
    [m,par0] = detrenddataNEW(raw,x,DetrendParam);
    y1 = y1-m; %detrended y1 
    y1 = y1 - mean(y1);
    

% fit OU and OUoscillatory models

    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise);
    par0M(i,:) = par0;
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
% show figure - comment out "showfigure" if don't want 
%individual cell figures popping up
    showfigure(x,m,raw,y1,BICdiff,par1,par2,i)

end
% updates total list for all pooled experiments
 BICdiffTOT1 = [BICdiffTOT1;BICdiffM];
 par0TOT1 = [par0TOT1;par0M];
 par1TOT1 = [par1TOT1;par1M];
 par2TOT1 = [par2TOT1;par2M];

end

BICdiffTOT = BICdiffTOT1;

%%
% plots distribution of LLR scores
figure()
hist(BICdiffTOT,15) % choose number of bins with second number
xlabel('LLR score')
ylabel('Frequency')
title('Distribution of LLR scores')
% xlim([-5,25])

%%
% calculate sample lengths to create synthetic data
sampTOT1 = [];
for j =1:length(fields)
j
current = data.(fields{j})
samp = [];
% load cell data for current experiment - loops through cells
for i = 1:size(current.data,2);
%     i
    y1 = current.data(:,i);  
    x = time;
    x(y1==0) = []; %deletes times from which no signal
    samp(i,1) = length(x)
    y1(y1==0) = [];
    y1 = y1 - mean(y1);
  
end
% updates total list for all pooled experiments
 sampTOT1 = [sampTOT1;samp];

end

repeats = round(2000/length(par1TOT1)); %to create roughly 2000 synthetic cells - this can be increased or decreased

%%
% turn on parallel - uncomment next line
% only need to run 'matlabpool open' once per session

% matlabpool open 4 % opens extra cores for parallel computing 

[ synthOUhier1 ] = MakesynthOUHIERACHICALvariablex(par0TOT1, par1TOT1,repeats,sampTOT1,time); % creates synthetic data from OU model

[ BICdiffsynthTOT ] = BICdistDATAnew( synthOUhier1,time,par1TOT1,repeats,DetrendParam); % find LLR distribution of snth OU data

%%
% LLRs can be tiny and just negative - this just sets them to zero
BICdiffTOT1(BICdiffTOT1<0) = 0;
BICdiffsynthTOT(BICdiffsynthTOT<0) = 0;

%%
% this sections calucaltes the FDR - first need to find pi_0
% pi_0 compares shape of data with the synth OU to estimate rough
% proportion of cells oscillating

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
    
% 
% figure()
% subplot(1,2,1)
% plot(range,piest)

% cubic spline regression
xx = linspace(lower,upper,100);
yy = spline(range,piest,xx);
% plot(xx,yy)

% %plot piGuess
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

%%
% Go through BIClist calculating q values

[BICdiffTOT,I] = sort(BICdiffTOT);

q1 = zeros(length(BICdiffTOT),1);
% find BIC list location above thresh...

for i = 1:length(BICdiffTOT)

    Thresh = BICdiffTOT(i);
    (sum(BICdiffsynthTOT>Thresh)/length(BICdiffsynthTOT));
    (sum(BICdiffTOT>Thresh)/length(BICdiffTOT));
    q1(i) = piGUESS1*(sum(BICdiffsynthTOT>Thresh)/length(BICdiffsynthTOT))/(sum(BICdiffTOT>Thresh)/length(BICdiffTOT));
end

q = 0.05; % choose initial cutoff level...
cutoff = find(q1<q,1,'first')
[w,l] = sort(I);
Reorderedq = q1(l);
PassList = Reorderedq<q;

%%
%only plot periods of those passing

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
hist(par2TOT(:,2)./par2TOT(:,1),[0:0.2:2])
xlim([0 3])
title('Quality of fitted oscillations to all cells')
xlabel('length scale (hours)')
ylabel('frequency')
%%
% plot LLR distributions of Hes1, control and bootstrap
subplot(1,2,1)
histogram(BICdiffTOT1(1:19),15)
t = title('Hes1 promoter');
t.FontWeight = 'normal';
xlim([0 41])
ylim([0 4])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


subplot(1,2,2)
histogram(BICdiffsynthTOT,15)
xlim([0 41])
t = title('Synthetic bootstrap (non-osc)');
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
