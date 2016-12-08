function [ synthOU ] = MakesynthOUHIERACHICALwithSE(par0, par1,repeats,x )
%SYNTHOU Summary of this function goes here
%   synthesises multiple OU traces

% create 100 OU and 100 OUosc data series - then pool

%% create 100 time series

cells = length(par0);
synthOU = [];

for i = 1:cells

x = x';
a0 = par0(i,1);%0.05; %0.1
b0 = par0(i,2);
a = par1(i,1);%0.05; %0.1
b = par1(i,2);
cov1 = b0*exp(-a0*x.^2) + b*exp(-a*x);
Noise = par1(i,3);
% cov2 = exp(-a*x).*cos(b*x);
% cov2 = exp(-0.05*x.^2);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


MU = zeros(repeats,length(x));
% Noise parameter represents standard dev of meas noise
% Noise = 0.1;
Meas = diag((Noise^2).*ones(1,length(x)));
CVM1 = CVM1 + Meas;
SIGMA = CVM1; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);

x = x';
% y1 = y1';
synthOUcurr = data1';
synthOU = [synthOU, synthOUcurr];
end

end

