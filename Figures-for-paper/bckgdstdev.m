function [stdev] = bckgdstdev(bckgd,time,lscale);
%BCKGDSTDEV Summary of this function goes here
%   Detailed explanation goes here

%%
% detrend with SE

for i = 1:size(bckgd,2)

y1 = bckgd(:,i);  
x = time;
x(y1==0) = []; %deletes times from which no signal
y1(y1==0) = [];
y1 = y1 - mean(y1);

likfunc = @likGauss; 
covfunc = @covSEa;
hyp2.cov = [lscale-1,log(std(y1))];
hyp2.lik = log(std(y1));

x1 = -10; x2 = lscale; eta = 16;
prior.cov = {{@priorSmoothBox2,x1,x2,eta}; {}};
inf = {@infPrior,@infExact,prior};
% inf = {@infExact,@infPrior,prior};

% hyp2 = minimize(hyp2, @gp, -10000, @infExact, [], covfunc, likfunc, x, y1);
hyp2 = minimize(hyp2, @gp, -10000, inf, [], covfunc, likfunc, x, y1);

nlmlOU = gp(hyp2, @infExact, [], covfunc, likfunc, x, y1);

stdvec(i) = exp(hyp2.lik);

clear prior

subplot(2,2,i)
z = x;%linspace(0, x(end), 1000)';
[m s2] = gp(hyp2, @infExact, [], covfunc, likfunc, x, y1, z);
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
fill([z; flipdim(z,1)], f, [7 7 7]/8);
hold on
plot(z,m)
plot(x,y1);%,'+')
hold off
xlabel('time (hours)')
ylabel('mean subtracted expression')
str = sprintf('St dev = %f',exp(hyp2.lik));
title(str);

end

stdev = mean(stdvec);

%%
end

