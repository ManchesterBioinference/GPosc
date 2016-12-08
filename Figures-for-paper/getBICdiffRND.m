function [BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise)
%GETBICDIFF Summary of this function goes here
%   loop through 10 times - find BIC with highest score
%%
avec = 0.001 + (0.5)*rand(10,1);
bvec = 2 + (4)*rand(10,1);
BICdiffvec = zeros(10,1);
par1vec = zeros(10,3);
par2vec = zeros(10,4);
%%
for i = 1:10;
    % create range to samplr from
    
samp = length(x);
likfunc = @likGauss; 
par1mean = [log(0.5),log(var(y1))];
par2mean = [log(avec(i)),log(2*pi/bvec(i)),log(var(y1))];
% par2mean = [log(0.001),log(2*pi/2),0]; only 6/1000 failed with these
% settings
covfunc = @covOUa; 
    hyp2.lik = log(Noise);
    hyp2.cov = par1mean;
    prior.lik = {{@priorDelta}};
    inf = {@infPrior,@infExact,prior};
    hyp2 = minimize(hyp2, @gp, -10000, inf, [], covfunc, likfunc, x, y1);
%     hyp2 = minimize(hyp2, @gp, -10000, @infExact, [], covfunc, likfunc, x, y1);
    nlmlOU = gp(hyp2, @infExact, [], covfunc, likfunc, x, y1);
    par1vec(i,:) = [exp(hyp2.cov),exp(hyp2.lik)];
    BIC1 = 2*nlmlOU/length(x);% + 2*log(samp);
    clear prior
    covfunc = @covOUosca; 
    hyp2.lik = log(Noise);
    hyp2.cov = par2mean;
    prior.lik = {{@priorDelta}};
    inf = {@infPrior,@infExact,prior};
    hyp2 = minimize(hyp2, @gp, -10000, inf, [], covfunc, likfunc, x, y1);
%     hyp2 = minimize(hyp2, @gp, -10000, @infExact, [], covfunc, likfunc, x, y1);
    nlmlOSC = gp(hyp2, @infExact, [], covfunc, likfunc, x, y1);
    par2vec(i,:) = [exp(hyp2.cov),exp(hyp2.lik)];
%     hyp2.lik
    BIC2 = 2*nlmlOSC/length(x);% + 3*log(samp);       
    BICdiffvec(i) = BIC1 - BIC2;
    clear prior
end
[BICdiff,posish] = max(BICdiffvec);
BICdiff = 100*BICdiff;
par1 = par1vec(posish,:);
par2 = par2vec(posish,:);
end

