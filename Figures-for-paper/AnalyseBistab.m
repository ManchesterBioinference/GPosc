par(1) = 2; 
par(2) = 0.03;
par(3) = 12;
options=ddeset('RelTol',1e-5);
tfin = 1500;
N = [1];% 50];
% Output_Times = linspace(5000.01,10000,4096);%5000:30:8000;%linspace(5000.01,15000,4096);%0:20:30000;
% Output_Times = linspace(0,tfin,tfin/30);%5000:30:8000;
Output_Times = 5000:30:(5000+tfin);
totalreps = 10;
for i=1:1
% check tobias to get approx values
mstart = 15*N;
pstart = 20*N;

tic
[ mout,pout] = GillespieBISTABmod(N,par,totalreps,mstart, pstart,Output_Times);

toc


end

% plot(Output_Times,mout1)
% subplot(2,3,2)
plot(Output_Times,mout,'color','b')
hold on
plot(Output_Times,pout,'color','r')
ylim([0 60*N])
xlabel('time (mins)')
ylabel('mRNA copy number')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'b)'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

%%
% fit OU and OUosc models
x = Output_Times';
% y1 = y1';
dataTOT = pout';
x = (x-5000)/60;
%%
% figure()
% plot(x,dataTOT(:,115))
% plot(x,cov2)

%%
Noise = sqrt(0.1);
dataNORMED = zeros(size(dataTOT));
% go through data and add noise
for i = 1:size(dataTOT,2);
%     i
    y1 = dataTOT(:,i);  
%     x(y1==0) = []; %deletes times from which no signal
    samp = length(x);
%     y1(y1==0) = [];
    y1 = y1 - mean(y1);
    
% remove trend from data

%     Noise = stdev/std(y1);
    y1 = y1/std(y1);
    MU = zeros(1,length(x));
%     Noise = 0.1;
    Meas = diag((Noise^2).*ones(1,samp));
    CVM1 = Meas;
    SIGMA = CVM1; % change this to switch non-osc and osc
    measerror = mvnrnd(MU,SIGMA);
    y2 = y1 + measerror';
    dataNORMED(:,i) = y2; 

end

%%


% fit models to data 

par1M = zeros(size(dataTOT,2),3);
par2M = zeros(size(dataTOT,2),4);
BICdiffM = zeros(size(dataTOT,2),1);

% load cell data for current experiment - loops through cells
for i = 1:size(dataTOT,2);
    i
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


%%

