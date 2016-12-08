%% 
% generate synthetic data
% verify with power spectrum that not oscillating
% 1000 non-osc and osc - generate from gillespie..

%%
%generate non-oscillatory time series
par(1) = 300;
par(2) = 1;
par(3) = 0.07;%0.03;
par(4) = 0.07;%0.03;
par(5) = 1;
par(6) = 1;
par(7) = 0;

%Approximate size used to preallocate space for particle number vector
% N = 50; %System size - sets total number of particles in system


Nvec = [20];
Output_Times = linspace(5000.01,10000,4096);
totalreps = 1000;
for i=1:1
mstart = 1*N;
pstart = 30*N;
tic
[ mout,pout] = GillespieTiming_Mod9nodelay(N,par,totalreps,mstart, pstart,Output_Times); 
toc
end

% figure()
% plot(Output_Times,mout(3,:))

for i = 1:totalreps
mfin=mout(i,:);
pfin=pout(i,:);%pout_temp=pout(i,:);    
fluc = (pfin/N - (mean(pfin))/N);
pf = (N^0.5)*fluc;
fluc = (mfin/N - (mean(mfin))/N);
mf = (N^0.5)*fluc;

Fs = 1/(Output_Times(2)-Output_Times(1));
L = length(pf);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y1 = fft(pf,NFFT)/L;
Y2 = fft(mf,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
f(1) = [];
% 
pre = (1/Fs*L)^0.5; 
FourierVecP1 = pre*(Y1(1:NFFT/2+1)); 
FourierVecM1 = pre*(Y2(1:NFFT/2+1));
PowerVecP(i,:) = FourierVecP1.*conj(FourierVecP1);
PowerVecM(i,:) = FourierVecM1.*conj(FourierVecM1);
CrossSpecMP(i,:) = FourierVecM1.*conj(FourierVecP1);
end

MeanPowerP = mean(PowerVecP);
MeanPowerM = mean(PowerVecM);
MeanCrossMP = mean(CrossSpecMP);

MeanPowerP(1) = [];
MeanPowerM(1) = [];
MeanCrossMP(1) = [];

fCAL = f;
MeanPowerOU = MeanPowerP;

%%
% generate oscillatory time series
par(1) = 100;
par(2) = 3;
par(3) = 0.03;
par(4) = 0.03;
par(5) = 1;
par(6) = 1;
par(7) = 18;

Nvec = [20];
Output_Times = linspace(5000.01,10000,4096);
totalreps = 1000;
for i=1:1
N = Nvec(i)
mstart = 1*N;
pstart = 30*N;

tic
[ mout,pout] = GillespieTiming_Mod9(N,par,totalreps,mstart, pstart,Output_Times); 
toc
end

% figure()
% plot(Output_Times,mout(3,:))

% confirm dynamical regime by plotting average power spec

for i = 1:totalreps
mfin=mout(i,:);
pfin=pout(i,:);  
fluc = (pfin/N - (mean(pfin))/N);
pf = (N^0.5)*fluc;
fluc = (mfin/N - (mean(mfin))/N);
mf = (N^0.5)*fluc;

Fs = 1/(Output_Times(2)-Output_Times(1));
L = length(pf);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y1 = fft(pf,NFFT)/L;
Y2 = fft(mf,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
f(1) = [];
% 
pre = (1/Fs*L)^0.5; 
FourierVecP1 = pre*(Y1(1:NFFT/2+1)); 
FourierVecM1 = pre*(Y2(1:NFFT/2+1));
PowerVecP(i,:) = FourierVecP1.*conj(FourierVecP1);
PowerVecM(i,:) = FourierVecM1.*conj(FourierVecM1);
CrossSpecMP(i,:) = FourierVecM1.*conj(FourierVecP1);
end

MeanPowerP = mean(PowerVecP);
MeanPowerM = mean(PowerVecM);
MeanCrossMP = mean(CrossSpecMP);

MeanPowerP(1) = [];
MeanPowerM(1) = [];
MeanCrossMP(1) = [];

MeanPowerOUosc = MeanPowerP;

%% 
% generate 100 time series for tests

%non osc
par(1) = 400;
par(2) = 5;
par(3) = 0.03;
par(4) = 0.03;
par(5) = 0.1;
par(6) = 1;
par(7) = 0;

Nvec = [20];
Output_Times = 5000:30:6500;
totalreps = 100;
for i=1:1
N = Nvec(i)
mstart = 1*N;
pstart = 30*N;

tic
[ mout1,pout1] = GillespieTiming_Mod9nodelay(N,par,totalreps,mstart, pstart,Output_Times); 

toc


end


par(1) = 100;
par(2) = 3;
par(3) = 0.03;
par(4) = 0.03;
par(5) = 1;
par(6) = 1;
par(7) = 18;

Nvec = [20];% 50];
for i=1:1

mstart = 1*N;
pstart = 30*N;

tic
[ mout2,pout2] = GillespieTiming_Mod9(N,par,totalreps,mstart, pstart,Output_Times); 

toc
end

%%
% fit OU and OUosc models
x = Output_Times';

data1 = pout1';
data2 = pout2';
dataTOT = [data1,data2];

likfunc = @likGauss; 

%%
x = x/60; % renormalise time to hours

%%
dataNORMED = zeros(size(dataTOT));
% go through data and add noise
for i = 1:size(dataTOT,2);
    i
    y1 = dataTOT(:,i);  
    samp = length(x);
    y1 = y1 - mean(y1);
    y1 = y1/std(y1);
    MU = zeros(1,length(x));
    Noise = 0.1;
    Meas = diag((Noise^2).*ones(1,samp));
    CVM1 = Meas;
    SIGMA = CVM1; % change this to switch non-osc and osc
    measerror = mvnrnd(MU,SIGMA);
    y2 = y1 + measerror';
    dataNORMED(:,i) = y2; 

end

%%
% % sanity check that adding noise worked..
% figure()
% plot(x,y1)
% figure()
% plot(x,y2)

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

%%
% plot(x,dataNORMED(:,1))
%%
% updates total list for all pooled experiments
 BICdiffTOT = BICdiffM;%[BICdiffTOT;BICdiffM];
 par1TOT = par1M;%[par1TOT;par1M];
 par2TOT = par2M;%[par2TOT;par2M];


%%

% load fig4

ratio = 1;%0.5/(sqrt(2*pi)); 

subplot(3,2,1)
plot(fCAL*60,log10(ratio*MeanPowerOU))
xlim([0 2])
xlabel('Frequency (1/hours)')
ylabel('log_{10}(Power)')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


subplot(3,2,2)
plot(fCAL*60,log10(ratio*MeanPowerOUosc))
xlim([0 2])
xlabel('Frequency (1/hours)')
ylabel('log_{10}(Power)')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

ex1 = 8; % 8 good
y1 = dataNORMED(:,ex1);

thr = 0.05;
% Pfa = [ 5 1 0.01]/100;
Pd = 1-thr;

[pxx,f,pth] = plomb(y1,x,'normalized','Pd',Pd);

sum(pxx>pth(1))
beat = (sum(pxx>pth(1))>0);

subplot(3,2,3)
plot(x-5000/60,y1)
xlabel('Time (hours)')
ylabel('Normalised protein')
number = 100*BICdiffM(ex1);
str = sprintf('LLR: %.2f',number);
title(str,'fontweight','normal');
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(3,2,5)
plot(f,pxx,f,pth*ones(size(f')))
xlabel('Frequency (1/hours)')
ylabel('Normalised power')
text(0.3*[1],pth-.5,[repmat('P_{fa} = ',[1 1]) num2str(thr')])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'E'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
timeseriesOU = y1;

ex = 110 %101 and 110 good exs;
y1 = dataNORMED(:,ex);

thr = 0.05;
% Pfa = [ 5 1 0.01]/100;
Pd = 1-thr;

[pxx,f,pth] = plomb(y1,x,'normalized','Pd',Pd);

sum(pxx>pth(1))
beat = (sum(pxx>pth(1))>0);

subplot(3,2,4)
plot(x-5000/60,y1)
xlabel('Time (hours)')
ylabel('Normalised protein')
number = 100*BICdiffM(ex);
str = sprintf('LLR: %.2f',number);
title(str,'fontweight','normal');
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'D'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(3,2,6)
plot(f,pxx,f,pth*ones(size(f')))
xlabel('Frequency (1/hours)')
ylabel('Normalised power')
text(0.3*[1],pth-.5,[repmat('P_{fa} = ',[1 1]) num2str(thr')])
timeseriesOUosc = y1;
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'F'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
