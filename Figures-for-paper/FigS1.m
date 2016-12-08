x = linspace(0,800,800);
%%
% define parameters for generating time series
% note all s.t.d kept at 1
par(1) = 0.15;
par(2) = 0.5;

%%
a = par(1);%0.05; %0.1
b = par(2);
Q2 = b/(2*pi()*a);
cov1 = exp(-a*x).*cos(b*x);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


MU = zeros(1,length(x));

SIGMA = CVM1; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);


OUosc = data1;
plot(x,data1)

%%
% load FigS1
plot(x,OUosc)
xlabel('Time (hours)')
ylabel('y(time)')
ylim([-4 4])