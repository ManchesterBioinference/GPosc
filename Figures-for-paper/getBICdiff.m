function [BICdiff, par1, par2] = getBICdiff(x,y1,Noise)
%GETBICDIFF Summary of this function goes here
%   Detailed explanation goes here
samp = length(x);
likfunc = @likGauss; 
par1mean = [log(0.5),0];
par2mean = [log(0.1),log(2*pi/6),0];
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
    par1 = [exp(hyp2.cov),exp(hyp2.lik)];
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
    par2 = [exp(hyp2.cov),exp(hyp2.lik)];
%     hyp2.lik
    BIC2 = 2*nlmlOSC/length(x);% + 3*log(samp);       
    BICdiff = BIC1 - BIC2;
    clear prior

end

