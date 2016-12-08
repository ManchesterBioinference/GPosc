function showLomb(x,m,raw,y1,i)
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
    title(str,'fontweight','normal');
    
    thr = 0.05;
    % Pfa = [ 5 1 0.01]/100;
    Pd = 1-thr;

    [pxx,f,pth] = plomb(y1,x,'normalized','Pd',Pd);

    sum(pxx>pth(1));
    beat = (sum(pxx>pth(1))>0);
    [M,I] = max(pxx);
    domfreq = 1/f(I);
    
    subplot(2,1,2)    
    plot(f,pxx,f,pth*ones(size(f')))
    xlabel('Frequency (1/hours)')
    ylabel('Normalised power')
    text(0.3*[1],pth-.5,[repmat('P_{fa} = ',[1 1]) num2str(thr')])
    a = xlim();
    b = ylim(); 
    if beat
    str = sprintf(['Dominant frequency = ',num2str(domfreq,2),' hours']);
    title(str,'fontweight','normal');
    else
%     str = sprintf('LLR score %f, Period = %f, Quality = %f',BICdiff,(2*pi()/par2(2)),par2(2)/par2(1));
%     title(str,'fontweight','normal');

end

