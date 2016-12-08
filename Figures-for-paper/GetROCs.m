function [ FP1, TP1, FP2, TP2 ] = getROCs( par1,par2, Tfinal, Noise, CellNum )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% 

Nvec = [20];% 50];
Output_Times = 5000:30:(5000+Tfinal);
totalreps = CellNum;
for i=1:1
N = Nvec(i)
% check tobias to get approx values
mstart = 1*N;
pstart = 30*N;

tic
[ mout1,pout1] = GillespieTiming_Mod9nodelay(N,par1,totalreps,mstart, pstart,Output_Times); 

toc


end


Nvec = [20];% 50];

for i=1:1
N = Nvec(i)
mstart = 1*N;
pstart = 30*N;

tic
[ mout2,pout2] = GillespieTiming_Mod9(N,par2,totalreps,mstart, pstart,Output_Times); 

toc


end

%%
% fit OU and OUosc models
x = Output_Times';
% y1 = y1';
data1 = pout1';
data2 = pout2';
dataTOT = [data1,data2];

likfunc = @likGauss; 

%%
x = x/60;

%%
dataNORMED = zeros(size(dataTOT));
% go through data and add noise
for i = 1:size(dataTOT,2);
    i
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

    [BICdiff, par1, par2] = getBICdiff(x,y1,Noise);
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;
% show figure - comment out "showfigure" is don't want 
%individual cell figures popping up
%     showfigure(x,m,raw,y1,BICdiff,par1,par2,i)

end

%%
% updates total list for all pooled experiments
 BICdiffTOT = BICdiffM;%[BICdiffTOT;BICdiffM];
 par1TOT = par1M;%[par1TOT;par1M];
 par2TOT = par2M;%[par2TOT;par2M];
 
 
 %%

% true +ve rate - false +ve rate - ROC curve

A =  BICdiffTOT(1:CellNum);
B = BICdiffTOT(CellNum+1:2*CellNum);


%%
% for given threshold - get sensitivity and specificity

thresh = linspace(min(BICdiffTOT)-1,max(BICdiffTOT)+1,200);

for i = 1:length(thresh)
    FP1(i) = sum(A>thresh(i));
    TP1(i) = sum(B>thresh(i));
    
end


%%

thrvec = exp(linspace(-15,-0.01,200));
beatvec = zeros(length(dataTOT),1);
%%

for i = 1:length(thrvec)
    thr = thrvec(i)
    for j = 1:length(dataTOT)
        y1 = dataTOT(:,j);
%         thr = 0.05;        
        Pd = 1-thr;
        [pxx,f,pth] = plomb(y1,x,'normalized','Pd',Pd);
        beatvec(j) = (sum(pxx>pth(1))>0);       
        
    end
    A =  beatvec(1:CellNum);
    B = beatvec(CellNum+1:2*CellNum);
    FP2(i) = sum(A);
    TP2(i) = sum(B);    
    
end


end

