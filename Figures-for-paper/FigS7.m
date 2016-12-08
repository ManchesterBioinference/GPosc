par1(1) = 300;
par1(2) = 1;
par1(3) = 0.07;%0.03;
par1(4) = 0.07;%0.03;
par1(5) = 1;
par1(6) = 1;
par1(7) = 0;


par2(1) = 100;
par2(2) = 4;%3;
par2(3) = 0.03;
par2(4) = 0.03;
par2(5) = 1;
par2(6) = 1;
par2(7) = 18;

%%

Tfinal = 1500; % 25 hours
Noise = sqrt(0.1); % was 0.1!!!!!
CellNum = 1000;

%Note have changed code - system size down!Default was 20
[ x, dataNORMED ] = GetTimeSeriesLow( par1,par2, Tfinal, Noise, CellNum );


%%


% fit models to data 

par1M = zeros(size(dataNORMED,2),3);
par2M = zeros(size(dataNORMED,2),4);
BICdiffM = zeros(size(dataNORMED,2),1);

% load cell data for current experiment - loops through cells
parfor i = 1:size(dataNORMED,2);
%     i
    y1 = dataNORMED(:,i);  
  

% fit OU and OUoscillatory models

    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise);
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;

end

% updates total list for all pooled experiments
 BICdiffTOT = BICdiffM;%[BICdiffTOT;BICdiffM];
 par1TOT = par1M;%[par1TOT;par1M];
 par2TOT = par2M;%[par2TOT;par2M];
 
 %%
repeats = 1; 
 % turn on parallel
[ synthOUhier1 ] = MakesynthOUHIERACHICAL( par1TOT,repeats,x );

[ BICdiffsynthTOT ] = BICdistDATAsynth( synthOUhier1,x,par1TOT,repeats);

save('FigS7')


%%

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
    
% cubic spline regression
xx = linspace(lower,upper,100);
yy = spline(range,piest,xx);
plot(xx,yy)

%plot piGuess
hold on
plot(xx,yy,'color','r')
ylim([0 1])
hold off
xlabel('\lambda')
ylabel('\pi_0(\lambda)')
xlim([lower1,upper])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'a)'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

piGUESS1 = min(1,yy(1));
%%

% Go through BIClist calculating q values
% choose initial cutoff level...

[BICdiffTOT,I] = sort(BICdiffTOT);

q1 = zeros(length(BICdiffTOT),1);
% find BIC list location above thresh...

for i = 1:length(BICdiffTOT)

    Thresh = BICdiffTOT(i);
    (sum(BICdiffsynthTOT>Thresh)/length(BICdiffsynthTOT));
    (sum(BICdiffTOT>=Thresh)/length(BICdiffTOT));
    q1(i) = piGUESS1*(sum(BICdiffsynthTOT>=Thresh)/length(BICdiffsynthTOT))/(sum(BICdiffTOT>=Thresh)/length(BICdiffTOT));
end

q = 0.05;%0.05;
cutoff = find(q1<q,1,'first')
[w,l] = sort(I);
Reorderedq = q1(l);
PassList = Reorderedq<q;
FP = sum(PassList(1:1000));
TP = sum(PassList(1001:2000));
FDR = FP/(FP+TP);

%%
cell = 1002;
y = dataNORMED;
subplot(1,2,1)
plot(x,y(:,cell))
xlabel('Time (hours)')
ylabel('Normalised expression')
a = xlim();
b = ylim();
text(a(1)-0.05*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(1,2,2)
histogram(y(:,cell))
xlim([-5 5])
xlabel('Normalised expression')
ylabel('Frequency')
a = xlim();
b = ylim();
text(a(1)-0.05*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

