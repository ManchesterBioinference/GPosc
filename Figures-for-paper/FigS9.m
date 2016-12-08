a1 = 0.001;
b1 = 2*pi()/24;
a2 = 0.1;
b2 = 2*pi()/2.5;
x = linspace(0,50,100)';
Noise = sqrt(0.1);

%% create 100 time series

cells = 10;
y = [];

for i = 1:cells

x = x';
% a1 = 1;%0.05; %0.1
% b1 = par0(i,2);
% a2 = par1(i,1);%0.05; %0.1
% b2 = par1(i,2);
cov1 = 5*exp(-a1*x).*cos(b1*x) + exp(-a2*x).*cos(b2*x);
% Noise = par1(i,3);
% cov2 = exp(-a*x).*cos(b*x);
% cov2 = exp(-0.05*x.^2);

CovMatrix1 = zeros(length(x),length(x));

for i = 1:length(x)
    CovMatrix1(:,i) = circshift(cov1',i-1);

end
    CovMatrix1 = CovMatrix1';
    
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';


MU = zeros(1,length(x));
% Noise parameter represents standard dev of meas noise
% Noise = 0.1;
Meas = diag((Noise^2).*ones(1,length(x)));
CVM1 = CVM1 + Meas;
SIGMA = CVM1; % change this to switch non-osc and osc
data1 = mvnrnd(MU,SIGMA);

x = x';
% y1 = y1';
synthOUcurr = data1';
y = [y, synthOUcurr];
end

%%
plot(x,y(:,3))
%%

% fit models to data 

par1M = zeros(size(y,2),3);
par2M = zeros(size(y,2),4);
BICdiffM = zeros(size(y,2),1);

% load cell data for current experiment - loops through cells
parfor i = 1:size(y,2);
    i
    y1 = y(:,i);  
  

% fit OU and OUoscillatory models

    [BICdiff, par1, par2] = getBICdiffRND24(x,y1,Noise);
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;

end

% updates total list for all pooled experiments
 BICdiffTOT1 = BICdiffM;%[BICdiffTOT;BICdiffM];
 par1TOT1 = par1M;%[par1TOT;par1M];
 par2TOT1 = par2M;%[par2TOT;par2M];
 
 %%
 
Lengthscale = 7.5; % roughly 3x period expected (in hours)
DetrendParam = log(1./(2*Lengthscale.^2));
 
par1M = zeros(size(y,2),3);
par2M = zeros(size(y,2),4);
BICdiffM = zeros(size(y,2),1);

% load cell data for current experiment - loops through cells
parfor i = 1:size(y,2);
    i
    y1 = y(:,i);  
    Noise1 = Noise/std(y1);
    y1 = y1/std(y1);
    raw = y1;
    [m,par0] = detrenddataNEW(raw,x,DetrendParam);
    y1 = y1-m; %detrended y1 
  

% fit OU and OUoscillatory models

    [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise1);
    par1M(i,:) = par1;
    par2M(i,:) = par2;
    BICdiffM(i,:) = BICdiff;

end

% updates total list for all pooled experiments
 BICdiffTOT2 = BICdiffM;%[BICdiffTOT;BICdiffM];
 par1TOT2 = par1M;%[par1TOT;par1M];
 par2TOT2 = par2M;%[par2TOT;par2M];
 
%%
y1 = y(:,1);  
Noise1 = Noise/std(y1);
y1 = y1/std(y1);
raw = y1;
[m,par0] = detrenddataNEW(raw,x,DetrendParam);
y1 = y1-m; %detrended y1 


% fit OU and OUoscillatory models

[BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise1);
par1M(i,:) = par1;
par2M(i,:) = par2;
BICdiffM(i,:) = BICdiff;
%%
subplot(2,1,1)
plot(x,raw)
xlabel('Time (hours)')
ylabel('Normalised expression')
hold on
plot(x,m,'r')
hold off
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(2,1,2)
plot(x,y1)
xlabel('Time (hours)')
ylabel('Normalised expression')
t = title(['Inferred period =',num2str(par2(2))]);
t.FontWeight = 'normal';
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.1*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')


 
 %%
repeats = 1; 
 % turn on parallel
[ synthOUhier1 ] = MakesynthOUHIERACHICAL( par1TOT1,repeats,x );

[ BICdiffsynthTOT1 ] = BICdistDATAsynth( synthOUhier1,x,par1TOT1,repeats);

[ synthOUhier2 ] = MakesynthOUHIERACHICAL( par1TOT2,repeats,x );

[ BICdiffsynthTOT2 ] = BICdistDATAsynth( synthOUhier2,x,par1TOT2,repeats);

