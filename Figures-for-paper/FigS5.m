% parameters for Gillespie

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
Noise1 = sqrt(0.1); % was 0.1!!!!!
CellNum = 1000;

[ x, dataNORMED ] = GetTimeSeries( par1,par2, Tfinal, Noise1, CellNum );


%%
y = dataNORMED;
currentdata = y;
par1M = zeros(size(currentdata,2),3);
par2M = zeros(size(currentdata,2),4);
BICdiffM = zeros(size(currentdata,2),1);
x = x';
time = x';

%%

% NO DETRENDING!

parfor i = 1:size(currentdata,2);
%     i
    y1 = currentdata(:,i);  
    x = time;
    x(y1==0) = []; %deletes times from which no signal
    samp = length(x);
    y1(y1==0) = [];
    y1 = y1 - mean(y1);
    
% remove trend from data

%     Noise = stdev/std(y1);
    Noise = Noise1/std(y1);
    y1 = y1/std(y1);
%     raw = y1;
% %IMPORTANT PARAMETER - number (3rd input) controls how slow
% % trend is for cell data
%     [m] = detrenddata(raw,x,-6);
%     y1 = y1-m; %detrended y1   
%     y1 = y1/std(y1);        

% fit OU and OUoscillatory models

    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise);
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
% show figure - comment out "showfigure" is don't want 
%individual cell figures popping up
%     showfigure(x,m,raw,y1,BICdiff,par1,par2,i)

end
 x = time;
 BICdiffTOT = [BICdiffM];
 par1TOT = [par1M];
 par2TOT = [par2M];
 

%%
% go through FDR calculation to check
% need parametric bootsrap for each
repeats = 1;
% turn on parallel
[ synthOUhier ] = MakesynthOUHIERACHICAL( par1TOT,repeats,x);
[ BICdiffsynthTOT ] = BICdistDATA( synthOUhier,time,par1TOT,repeats);

% 
% h = kstest2(BICdiffTOT6(1:1000),BICdiffsynthTOT6(1:1000))

disp('done Gillespienotrend')
% save('Gillespienotrend.mat') 
clear

%%


load FigS5

subplot(1,3,1)
% subplot(3,1,1)
hist(BICdiffTOT(1:1000),0.5:30)
t = title(['Gillespie. Mean = ',num2str(mean(BICdiffTOT(1:1000)),'%.2f')]);
t.FontWeight = 'normal';
xlim([0 30])
xlabel('LLR')
ylabel('Frequency')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(1,3,2)
% subplot(3,1,2)
hist(BICdiffsynthTOT(1:1000),0.5:30)
xlim([0 30])
xlabel('LLR')
ylabel('Frequency')
t = title(['Bootstrap. Mean = ',num2str(mean(BICdiffsynthTOT(1:1000)),'%.2f')]);
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(1,3,3)
% subplot(3,1,3)
qqplot((BICdiffsynthTOT(1:1000)),(BICdiffTOT(1:1000)))
xlabel('Quantiles of bootstrap OU data')
ylabel('Quantiles of Gillespie LLRs')
xlim([0 34])
ylim([0 34])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
