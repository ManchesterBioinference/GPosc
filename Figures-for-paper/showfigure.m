function showfigure(x,m,raw,y1,BICdiff,par1,par2,i)
%SHOWFIGURE Summary of this function goes here
%   Detailed explanation goes here

%%
    % plots raw data with trend and detrended data
    figure()
    subplot(2,1,1)
%     f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
%     fill([z; flipdim(z,1)], f, [7 7 7]/8);
    hold on
    plot(x,m,'color','r')
    plot(x,raw);%,'+')
    hold off
    xlabel('Time (hours)')
    str = sprintf('Raw expression cell %.0f',i);
    title(str);
    subplot(2,1,2)    
    plot(x,y1)
    xlabel('Time (hours)')
    if BICdiff<3
    str = sprintf('BIC score %f',BICdiff);
    title(str);
    else
    str = sprintf('BIC score %f, Period = %f, Quality = %f',BICdiff,(2*pi()/par2(2)),1/par2(1));
    title(str);        

end

