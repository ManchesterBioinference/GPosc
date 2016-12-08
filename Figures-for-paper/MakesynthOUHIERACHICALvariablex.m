function [ synthOU ] = MakesynthOUHIERACHICALvariablex( par0, par1,repeats,sampTOT,time )
%SYNTHOU Summary of this function goes here
%   synthesises multiple OU traces

%create empty matrix to store outputs in
% XVec=repmat(sampTOT,1,repeats)';XVec=XVec(:);

% create 100 OU and 100 OUosc data series - then pool
% x = (0:0.5:avelength/2)';
%% create 100 time series

cells = length(par1);
synthOU = [];

for i = 1:cells

% x = x';
x = time(1:sampTOT(i))';
% disp(x)

a0 = par0(i,1);%0.05; %0.1
b0 = par0(i,2);
a = par1(i,1);%0.05; %0.1
b = par1(i,2);
cov1 = b0*exp(-a0*x.^2) + b*exp(-a*x);
Noise = par1(i,3);

CovMatrix1 = zeros(length(x),length(x));

for j = 1:length(x)
    CovMatrix1(:,j) = circshift(cov1',j-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


MU = zeros(repeats,length(x));
% Noise parameter represents standard dev of meas noise
% Noise = 0.1;
Meas = diag((Noise^2).*ones(1,length(x)));
CVM1 = CVM1 + Meas;
SIGMA = CVM1; % change this to switch non-osc and osc
% disp(MU)
% disp(SIGMA)
data1 = mvnrnd(MU,SIGMA);

% x = x';
% y1 = y1';
% synthOUcurr = data1';

synthOUtemp = zeros(max(sampTOT),repeats);
synthOUtemp(1:size(data1,2),:) = data1';
% disp(synthOUtemp)
synthOU = [synthOU, synthOUtemp];
end

end