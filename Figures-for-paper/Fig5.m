% run startup

par1(1) = 300;
par1(2) = 1;
par1(3) = 0.07;%0.03;
par1(4) = 0.07;%0.03;
par1(5) = 1;
par1(6) = 1;
par1(7) = 0;

par2(1) = 100;
par2(2) = 3;
par2(3) = 0.03;
par2(4) = 0.03;
par2(5) = 1;
par2(6) = 1;
par2(7) = 18;

%%

Tfinal = 1500; % 25 hours
Noise = sqrt(0.1);
CellNum = 1000;


[ FP11, TP11, FP21, TP21 ] = GetROCs( par1,par2, Tfinal, Noise, CellNum )


%%
Noise = sqrt(0.5);

[ FP12, TP12, FP22, TP22 ] = GetROCs( par1,par2, Tfinal, Noise, CellNum )

%%
Noise = sqrt(0.1);
Tfinal = 600; %10 hours

[ FP13, TP13, FP23, TP23 ] = GetROCs( par1,par2, Tfinal, Noise, CellNum )

% load fig5
%%
figure()
plot(FP11/CellNum,TP11/CellNum,'color','r')
hold on
plot(FP21/CellNum,TP21/CellNum,'color','b')
xlabel('1 - Specificity (false positive rate)')
ylabel('Sensitivity (true positive rate)')


%%
plot(FP12/CellNum,TP12/CellNum,'color','r')
hold on
plot(FP22/CellNum,TP22/CellNum,'color','b')
xlabel('1 - Specificity (false positive rate)')
ylabel('Sensitivity (true positive rate)')

%%

figure()
plot(FP13/CellNum,TP13/CellNum,'color','r')
hold on
plot(FP23/CellNum,TP23/CellNum,'color','b')
xlabel('1 - Specificity (false positive rate)')
ylabel('Sensitivity (true positive rate)')

%%

% create the associated time series
[ x1, dataNORMED1 ] = GetTimeSeries( par1,par2, 1500, sqrt(0.1) );

% create the associated time series
[ x2, dataNORMED2 ] = GetTimeSeries( par1,par2, 1500, sqrt(0.5) );

[ x3, dataNORMED3 ] = GetTimeSeries( par1,par2, 600, sqrt(0.1) );


%%

subplot(3,3,1)
plot(x1,dataNORMED1(:,1))
xlim([ 0 max(x1)])
xlabel('Time (hours)')
ylabel('Normalised expression')
ylim([-4 4])
text(-15, 1, '\sigma_N^{2} = 0.1')
text(-15, -1, '25 hours')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'A'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(3,3,2)
plot(x1,dataNORMED1(:,2))
xlim([ 0 max(x1)])
xlabel('Time (hours)')
ylabel('Normalised expression')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'B'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
subplot(3,3,3)
plot(FP11/CellNum,TP11/CellNum,'color','r')
hold on
plot(FP21/CellNum,TP21/CellNum,'color','b')
xlabel('1 - Specificity (false positive rate)')
ylabel('Sensitivity (true positive rate)')
legend('GP','L-S','Location','southeast')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'C'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(3,3,4)
plot(x2,dataNORMED2(:,1))
xlim([ 0 max(x2)])
xlabel('Time (hours)')
ylabel('Normalised expression')
ylim([-4 4])
text(-15, 1, '\sigma_N^{2} = 0.5')
text(-15, -1, '25 hours')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'D'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(3,3,5)
plot(x2,dataNORMED2(:,2))
xlim([ 0 max(x2)])
xlabel('Time (hours)')
ylabel('Normalised expression')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'E'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(3,3,6)
plot(FP12/CellNum,TP12/CellNum,'color','r')
hold on
plot(FP22/CellNum,TP22/CellNum,'color','b')
xlabel('1 - Specificity (false positive rate)')
ylabel('Sensitivity (true positive rate)')
legend('GP','L-S','Location','southeast')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'F'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(3,3,7)
plot(x3,dataNORMED3(:,1))
xlim([ 0 max(x3)])
xlabel('Time (hours)')
ylabel('Normalised expression')
ylim([-4 4])
text(-15/2.5, 1, '\sigma_N^{2} = 0.1')
text(-15/2.5, -1, '10 hours')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'G'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(3,3,8)
plot(x3,dataNORMED3(:,2))
xlim([ 0 max(x3)])
xlabel('Time (hours)')
ylabel('Normalised expression')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'H'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')

subplot(3,3,9)
plot(FP13/CellNum,TP13/CellNum,'color','r')
hold on
plot(FP23/CellNum,TP23/CellNum,'color','b')
xlabel('1 - Specificity (false positive rate)')
ylabel('Sensitivity (true positive rate)')
legend('GP','L-S','Location','southeast')
a = xlim();
b = ylim();
text(a(1)-0.1*(a(2)-a(1)),b(2)+0.03*(b(2)-b(1)),{'I'},...
    'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
