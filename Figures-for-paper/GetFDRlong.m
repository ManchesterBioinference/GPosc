function [ FP, TP, FDR ] = GetFDRlong( BICdiffTOT, BICdiffsynthTOT)
%GETFDR Summary of this function goes here
%   Detailed explanation goes here

upper = max([BICdiffTOT;BICdiffsynthTOT]);
lower1 = min([BICdiffTOT;BICdiffsynthTOT]);
lower = upper - 0.9*(upper-lower1);

range = linspace(lower,upper,20);

for i= 1:length(range)
    cutoff = range(i);
    num = sum(BICdiffTOT<cutoff)/length(BICdiffTOT);
    denom = sum(BICdiffsynthTOT<cutoff)/length(BICdiffsynthTOT);
    piest(i) =  num/denom;
end

%%
% cubic spline regression
xx = linspace(lower,upper,100);
yy = spline(range,piest,xx);
% plot(xx,yy)

% p = polyfit(range,piest,3);
% yy = polyval(p,xx);
    
%%
% figure()
% subplot(1,2,1)
% plot(range,piest)
% hold on
% plot(xx,yy,'color','r')
% ylim([0 1])
% hold off
% xlabel('\lambda')
% ylabel('\pi_0(\lambda)')
% xlim([lower1,upper])
% a = xlim();
% b = ylim();
% text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'a)'},...
%     'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
%%
piGUESS1 = yy(1);


% Go through BIClist calculating q values
% choose initial cutoff level...

[BICdiffTOT,I] = sort(BICdiffTOT);
% BICdiffsynthTOT = BICdiffsynthTOT6;


q1 = zeros(length(BICdiffTOT),1);
% find BIC list location above thresh...

for i = 1:length(BICdiffTOT)

    Thresh = BICdiffTOT(i);
    (sum(BICdiffsynthTOT>Thresh)/length(BICdiffsynthTOT));
    (sum(BICdiffTOT>Thresh)/length(BICdiffTOT));
    q1(i) = piGUESS1*(sum(BICdiffsynthTOT>Thresh)/length(BICdiffsynthTOT))/(sum(BICdiffTOT>Thresh)/length(BICdiffTOT));
end

q = 0.05;%0.05;
cutoff = find(q1<q,1,'first');
[~,l] = sort(I);
Reorderedq = q1(l);
PassList = Reorderedq<q;
FP = sum(PassList(1:2000))/length(PassList(1:2000));
TP = sum(PassList(2001:4000))/length(PassList(2001:4000));
FDR = FP/(FP+TP);

end
