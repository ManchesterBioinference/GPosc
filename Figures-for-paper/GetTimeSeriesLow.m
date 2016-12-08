function [ x, dataNORMED ] = GetTimeSeriesLow( par1,par2, Tfinal, Noise, totalreps )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Nvec = 1;%[20];% 50];
% Output_Times = linspace(5000.01,10000,4096);%5000:30:8000;%linspace(5000.01,15000,4096);%0:20:30000;
Output_Times = 5000:30:(5000+Tfinal);
% totalreps = 1;
for i=1:1
N = Nvec(i);
% check tobias to get approx values
mstart = 1*N;
pstart = 30*N;

% tic
[ mout1,pout1] = GillespieTiming_Mod9nodelay(N,par1,totalreps,mstart, pstart,Output_Times); 

% toc


end

for i=1:1
N = Nvec(i);
% check tobias to get approx values
mstart = 1*N;
pstart = 30*N;

% tic
[ mout2,pout2] = GillespieTiming_Mod9(N,par2,totalreps,mstart, pstart,Output_Times); 

% toc


end

%%
% fit OU and OUosc models
x = Output_Times';
% y1 = y1';
data1 = pout1';
data2 = pout2';
dataTOT = [data1,data2];


%%
x = (x-5000)/60;
%%
% figure()
% plot(x,dataTOT(:,115))
% plot(x,cov2)

%%
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

end

