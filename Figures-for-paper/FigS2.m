load Fig6

BICdiffTOT = (BICdiffTOT6);
BICdiffsynthTOT = (BICdiffsynthTOT6);

BICdiffTOT(BICdiffTOT<0) = 0;
BICdiffsynthTOT(BICdiffsynthTOT<0) = 0;

upper = max([BICdiffTOT;BICdiffsynthTOT]);
lower1 = min([BICdiffTOT;BICdiffsynthTOT]);
lower = upper - 0.9*(upper-lower1);



range = linspace(lower,upper,20);

for i= 1:length(range)
%     i
    cutoff = range(i);
    num = sum(BICdiffTOT<cutoff)/length(BICdiffTOT);
    denom = sum(BICdiffsynthTOT<cutoff)/length(BICdiffsynthTOT);
    piest(i) =  num/denom;
end
    

figure()
subplot(1,3,1)
scatter(range,piest,'x')

% cubic spline regression
xx = linspace(lower,upper,100);
yy = spline(range,piest,xx);
% plot(xx,yy)
% plot(range, piest)
% figure()

hold on
plot(xx,yy,'color','r')
ylim([0 1])
hold off
xlabel('\lambda')
ylabel('\pi_0(\lambda)')
xlim([lower1,upper])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
title('Trend exp(-6)','fontweight','normal')

piGUESS1 = yy(1);




BICdiffTOT = (BICdiffTOT4);
BICdiffsynthTOT = (BICdiffsynthTOT4);

BICdiffTOT(BICdiffTOT<0) = 0;
BICdiffsynthTOT(BICdiffsynthTOT<0) = 0;

upper = max([BICdiffTOT;BICdiffsynthTOT]);
lower1 = min([BICdiffTOT;BICdiffsynthTOT]);
lower = upper - 0.9*(upper-lower1);



range = linspace(lower,upper,20);

for i= 1:length(range)
%     i
    cutoff = range(i);
    num = sum(BICdiffTOT<cutoff)/length(BICdiffTOT);
    denom = sum(BICdiffsynthTOT<cutoff)/length(BICdiffsynthTOT);
    piest(i) =  num/denom;
end
    

% figure()
subplot(1,3,2)
scatter(range,piest,'x')

% cubic spline regression
xx = linspace(lower,upper,100);
yy = spline(range,piest,xx);
% plot(xx,yy)
% plot(range, piest)
% figure()

hold on
plot(xx,yy,'color','r')
ylim([0 1])
hold off
xlabel('\lambda')
ylabel('\pi_0(\lambda)')
xlim([lower1,upper])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
title('Trend exp(-4)','fontweight','normal')

piGUESS1 = yy(1);




BICdiffTOT = (BICdiffTOT2);
BICdiffsynthTOT = (BICdiffsynthTOT2);

BICdiffTOT(BICdiffTOT<0) = 0;
BICdiffsynthTOT(BICdiffsynthTOT<0) = 0;

upper = max([BICdiffTOT;BICdiffsynthTOT]);
lower1 = min([BICdiffTOT;BICdiffsynthTOT]);
lower = upper - 0.9*(upper-lower1);



range = linspace(lower,upper,20);

for i= 1:length(range)
%     i
    cutoff = range(i);
    num = sum(BICdiffTOT<cutoff)/length(BICdiffTOT);
    denom = sum(BICdiffsynthTOT<cutoff)/length(BICdiffsynthTOT);
    piest(i) =  num/denom;
end
    

% figure()
subplot(1,3,3)
scatter(range,piest,'x')

% cubic spline regression
xx = linspace(lower,upper,100);
yy = spline(range,piest,xx);
% plot(xx,yy)
% plot(range, piest)
% figure()

hold on
plot(xx,yy,'color','r')
ylim([0 1])
hold off
xlabel('\lambda')
ylabel('\pi_0(\lambda)')
xlim([lower1,upper])
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
title('Trend exp(-2)','fontweight','normal')

piGUESS1 = yy(1);


