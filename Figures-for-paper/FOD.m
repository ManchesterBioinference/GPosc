Lengthscale = 7.5; % roughly 3x period expected (in hours)
DetrendParam = log(1./(2*Lengthscale.^2));

%%

% run startup
% load in data from Excel file
num = xlsread('Hes1ublucF5_sparse_combine_30minexp_forBayesian.xlsx',1);
num(isnan(num)) = 0;


% assigns time 
% NOTE - important that time is in hours!!! Check
time = num(:,1);
time = time/(60*60*1000); %convert from ms to hours
time1 = time; % store for later...



% Manually enter the columns of data and associated bckgd for each experiment
data.exp1.data = num(:,2:8); data.exp1.bckgd = num(:,9:12);

data.exp2.data = num(:,13:24); data.exp2.bckgd = num(:,25:28);

% data.exp3.data = num(:,29:55); data.exp3.bckgd = num(:,56:59);

% data.exp4.data = num(:,23); data.exp4.bckgd = num(:,25:28);
% 
% data.exp5.data = num(:,23); data.exp5.bckgd = num(:,25:28);


%%
% creates empty vectors where final results will be stored
BICdiffTOT1 = [];
par0TOT1 = [];
par1TOT1 = [];
par2TOT1 = [];

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

par0M = zeros(size(current.data,2),3);
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
% show figure - comment out "showfigure" is don't want 
%individual cell figures popping up
%     showfigure(x,m,raw,y1,BICdiff,par1,par2,i)

end
% updates total list for all pooled experiments
 BICdiffTOT1 = [BICdiffTOT1;BICdiffM];
 par0TOT1 = [par0TOT1;par0M];
 par1TOT1 = [par1TOT1;par1M];
 par2TOT1 = [par2TOT1;par2M];

end

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


%%
% insert control data section
% laod control data
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
par0TOT2 = [];
par1TOT2 = [];
par2TOT2 = [];

% gets number of experiments
fields = fieldnames(data)

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

[stdev] = bckgdstdev(bckgd,time,-7)

par0M = zeros(size(current.data,2),3);
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
    [m,par0] = detrenddataNEW(raw,x,DetrendParam);
    y1 = y1-m; %detrended y1   
    y1 = y1 - mean(y1);
     
% fit OU and OUoscillatory models

    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise);
    par0M(i,:) = par0;
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
% show figure - comment out "showfigure" is don't want 
%individual cell figures popping up
%     showfigure(x,m,raw,y1,BICdiff,par1,par2,i)

end
% updates total list for all pooled experiments
 BICdiffTOT2 = [BICdiffTOT2;BICdiffM];
 par0TOT2 = [par0TOT2;par0M];
 par1TOT2 = [par1TOT2;par1M];
 par2TOT2 = [par2TOT2;par2M];

end


%%
% calculate sample lengths for control data...

sampTOT2 = [];
for j =1:length(fields)
% j
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
 sampTOT2 = [sampTOT2;samp];

end


%%
% add together 
BICdiffTOT1 = [BICdiffTOT1;BICdiffTOT2];
par0TOT1 = [par0TOT1;par0TOT2];
par1TOT1 = [par1TOT1;par1TOT2];
par2TOT1 = [par2TOT1;par2TOT2];
sampTOT1 = [sampTOT1;sampTOT2];
time = time1; %reset time to C17 data


repeats = round(2000/length(par1TOT1)); %to create roughly 2000 synthetic cells

%%
% turn on parallel
[ synthOUhier1 ] = MakesynthOUHIERACHICALvariablex(par0TOT1, par1TOT1,repeats,sampTOT1,time);

[ BICdiffsynthTOT ] = BICdistDATAnew( synthOUhier1,time,par1TOT1,repeats,DetrendParam);

%%
BICdiffTOT1(BICdiffTOT1<0) = 0;
BICdiffsynthTOT(BICdiffsynthTOT<0) = 0;


% BICdiffTOT = (BICdiffTOT1(20:44)); % for control data
BICdiffTOT = (BICdiffTOT1(1:19)); % for Hes1 data
BICdiffsynthTOT = (BICdiffsynthTOT);

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

q = 0.05;%0.05;
cutoff = find(q1<q,1,'first')
[w,l] = sort(I);
Reorderedq = q1(l);
PassList = Reorderedq<q;

%%
% plot LLR distributions of Hes1, control and bootstrap
subplot(1,3,1)
histogram(BICdiffTOT1(1:19),0:1:42)
t = title('Hes1 promoter');
t.FontWeight = 'normal';
xlim([0 41])
ylim([0 4])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(1,3,2)
histogram(BICdiffTOT1(20:44),0:1:41)
xlim([0 41])
ylim([0 11])
t = title('MMLV promoter');
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(1,3,3)
histogram(BICdiffsynthTOT,0:41)
xlim([0 41])
t = title('Synthetic bootstrap (non-osc)');
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


%%
% 
% 
% BICdiffTOT1(BICdiffTOT1<0) = 0;
% 
% BicTab = [];
% for i = 1:9
%     BicTab = [BicTab;'{\color{red} cell ' num2str(i) ' } '];
% end
% 
% for i = 10:19
%     BicTab = [BicTab;'{\color{red} cell ' num2str(i) '} '];
% end
% 
% for i = 20:44
%     BicTab = [BicTab;'{\color{blue} cell ' num2str(i) '}'];
% end
% 
% C = cellstr(BicTab);
% %%
% 
% [d1,d2] = sort(BICdiffTOT1(:,1),'descend');
% % d2 gives indexes of sorted BIC list
% BIClist = BicTab(d2,:); 
% BIClist = cellstr(BIClist);
% BICsorted = BICdiffTOT1(d2,:); 
% % C2 = num2cell(BICsorted);
% C2 = num2str(BICsorted);
% C2 = cellstr(C2);
% C1 = {C2,C};
% % sprintf('{\color{red} cell %d }', 1:19);
% 
% % BicTab = ['{\color{red} cell ' num2str(1) '}']
% 
% % make latex file
% rowNames = {'BIC score'; 'Cell number'};
% var1= C2;
% var2= BIClist;
% T = table(var1, var2);%,  'ColumnNames', rowNames);
% 
% % Now use this table as input in our input struct:
% input.data = T;
% 
% % Switch to generate a complete LaTex document or just a table:
% input.makeCompleteLatexDocument = 1;
% 
% % Switch transposing/pivoting your table if needed:
% input.transposeTable = 0;
% 
% % Switch to generate a complete LaTex document or just a table:
% input.makeCompleteLatexDocument = 1;
% 
% % Now call the function to generate LaTex code:
% latex = latexTable(input)

%%
% Make table of LLR values with -4.5, -5 and -4

BicTab = [];
for i = 1:9
    BicTab = [BicTab;'{\color{red} cell ' num2str(i) ' } '];
end

for i = 10:19
    BicTab = [BicTab;'{\color{red} cell ' num2str(i) '} '];
end

for i = 20:44
    BicTab = [BicTab;'{\color{blue} cell ' num2str(i) '}'];
end

C = cellstr(BicTab);
%%
load data45
BICdiffTOT1(BICdiffTOT1<0) = 0;
[d1,d2] = sort(BICdiffTOT1(:,1),'descend');
% d2 gives indexes of sorted BIC list
BIClist = BicTab(d2,:); 
BIClist = cellstr(BIClist);
BICsorted = BICdiffTOT1(d2,:); 
% C2 = num2cell(BICsorted);
C2 = num2str(BICsorted);
C2_1 = cellstr(C2);
% C1 = {C2,C};
% C_1 = C1;
BIClist_1 = BIClist;

load data4
BICdiffTOT1(BICdiffTOT1<0) = 0;
[d1,d2] = sort(BICdiffTOT1(:,1),'descend');
% d2 gives indexes of sorted BIC list
BIClist = BicTab(d2,:); 
BIClist = cellstr(BIClist);
BICsorted = BICdiffTOT1(d2,:); 
% C2 = num2cell(BICsorted);
C2 = num2str(BICsorted);
C2_2 = cellstr(C2);
% C1 = {C2,C};
% C_2 = C1;
BIClist_2 = BIClist;

load data5
BICdiffTOT1(BICdiffTOT1<0) = 0;
[d1,d2] = sort(BICdiffTOT1(:,1),'descend');
% d2 gives indexes of sorted BIC list
BIClist = BicTab(d2,:); 
BIClist = cellstr(BIClist);
BICsorted = BICdiffTOT1(d2,:); 
% C2 = num2cell(BICsorted);
C2 = num2str(BICsorted);
C2_3 = cellstr(C2);
% C1 = {C2,C};
% C_3 = C1;
BIClist_3 = BIClist;
% sprintf('{\color{red} cell %d }', 1:19);

% BicTab = ['{\color{red} cell ' num2str(1) '}']

% make latex file
rowNames = {'LLR score'; 'Cell number';'LLR score'; 'Cell number';'LLR score'; 'Cell number'};
var1= C2_1;
var2= BIClist_1;
var3= C2_2;
var4= BIClist_2;
var5= C2_3;
var6= BIClist_3;
T = table(var1, var2, var3, var4, var5, var6);%,  'ColumnNames', rowNames);

%%

% Now use this table as input in our input struct:
input.data = T;

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;

% Switch transposing/pivoting your table if needed:
input.transposeTable = 0;

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;

% Now call the function to generate LaTex code:
latex = latexTable(input)

