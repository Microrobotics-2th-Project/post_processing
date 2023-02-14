%Error group analysis
clear all
close all
myfig

Error_x_big_noflux = [];
Error_y_big_noflux = [];
Error_x_big_flux = [];
Error_y_big_flux = [];
Error_x_small_noflux = [];
Error_y_small_noflux = [];
Error_x_small_flux = [];
Error_y_small_flux = [];

%% load experiment error data
load('Error_exp1_bignoflux_30_06.mat')
Error_x_big_noflux = [Error_x_big_noflux; error_x];
Error_y_big_noflux = [Error_y_big_noflux, error_y];
load('Error_exp2_bignoflux_30_06.mat')
Error_x_big_noflux = [Error_x_big_noflux; error_x];
Error_y_big_noflux = [Error_y_big_noflux, error_y];
load('Error_exp3_bigflux_30_06.mat')
Error_x_big_flux = [Error_x_big_flux; error_x];
Error_y_big_flux = [Error_y_big_flux, error_y];
load('Error_exp4_smallflux_30_06.mat')
Error_x_small_flux = [Error_x_small_flux; error_x];
Error_y_small_flux = [Error_y_small_flux, error_y];
load('Error_exp5_smallnoflux_30_06.mat')
Error_x_small_noflux = [Error_x_small_noflux; error_x];
Error_y_small_noflux = [Error_y_small_noflux, error_y];
load('Error_exp6_smallnoflux_30_06.mat')
Error_x_small_noflux = [Error_x_small_noflux; error_x];
Error_y_small_noflux = [Error_y_small_noflux, error_y];

%% scatter plot of groups
figure
plot(Error_x_big_noflux,Error_y_big_noflux,'bo')
hold on
grid on
plot(Error_x_big_flux,Error_y_big_flux,'b^')
plot(Error_x_small_noflux,Error_y_small_noflux,'rs')
plot(Error_x_small_flux,Error_y_small_flux,'r*')
xlabel('Error x [mm]')
ylabel('Error y [mm]')
legend('Big MR no flux','Big MR flux','Small MR no flux','Small MR flux')

figure
plot(Error_x_big_noflux,Error_y_big_noflux,'bo')
hold on
grid on
plot(Error_x_small_noflux,Error_y_small_noflux,'ro')
xlabel('Error x [mm]')
ylabel('Error y [mm]')
legend('Big MR no flux','Small MR no flux')

figure
MR_diam_big = 1.25;
MR_diam_small = 0.55;
plot(Error_x_big_noflux./MR_diam_big,Error_y_big_noflux./MR_diam_big,'bo')
hold on
grid on
plot(Error_x_big_flux./MR_diam_big,Error_y_big_flux./MR_diam_big,'b^')
plot(Error_x_small_noflux./MR_diam_small,Error_y_small_noflux./MR_diam_small,'ro')
plot(Error_x_small_flux./MR_diam_small,Error_y_small_flux./MR_diam_small,'r^')
xlabel('Normalized error x []')
ylabel('Normalized error y []')
legend('Big MR no flux','Big MR flux','Small MR no flux','Small MR flux')

%% BoxPlot of groups
T = NaN(110,4);
% T(1:length(Error_x_big_noflux),1) = Error_x_big_noflux./MR_diam_big;
% T(1:length(Error_x_big_flux),2) = Error_x_big_flux./MR_diam_big;
% T(1:length(Error_x_small_noflux),3) = Error_x_small_noflux./MR_diam_small;
% T(1:length(Error_x_small_flux),4) = Error_x_small_flux./MR_diam_small;
T(1:length(Error_x_big_noflux),1) = Error_x_big_noflux;
T(1:length(Error_x_big_flux),2) = Error_x_big_flux;
T(1:length(Error_x_small_noflux),3) = Error_x_small_noflux;
T(1:length(Error_x_small_flux),4) = Error_x_small_flux;
figure
subplot(1,2,1)
boxplot(T,'Labels',{'Big no flux','Big flux','Small no flux','Small flux'})
grid on
ylabel('[mm]')
ylim([0 4])
title('Error on x')
T = NaN(110,4);
% T(1:length(Error_y_big_noflux),1) = Error_y_big_noflux./MR_diam_big;
% T(1:length(Error_y_big_flux),2) = Error_y_big_flux./MR_diam_big;
% T(1:length(Error_y_small_noflux),3) = Error_y_small_noflux./MR_diam_small;
% T(1:length(Error_y_small_flux),4) = Error_y_small_flux./MR_diam_small;
T(1:length(Error_y_big_noflux),1) = Error_y_big_noflux;
T(1:length(Error_y_big_flux),2) = Error_y_big_flux;
T(1:length(Error_y_small_noflux),3) = Error_y_small_noflux;
T(1:length(Error_y_small_flux),4) = Error_y_small_flux;
subplot(1,2,2)
boxplot(T,'Labels',{'Big no flux','Big flux','Small no flux','Small flux'})
grid on
ylabel('[mm]')
ylim([0 4])
title('Error on y')


T = NaN(110,4);
% T(1:length(Error_x_big_noflux),1) = Error_x_big_noflux./MR_diam_big;
% T(1:length(Error_x_big_flux),2) = Error_x_big_flux./MR_diam_big;
% T(1:length(Error_x_small_noflux),3) = Error_x_small_noflux./MR_diam_small;
% T(1:length(Error_x_small_flux),4) = Error_x_small_flux./MR_diam_small;
T(1:length(Error_x_big_noflux),1) = Error_x_big_noflux;
T(1:length(Error_x_big_flux),2) = -10;
T(1:length(Error_x_small_noflux),3) = Error_x_small_noflux;
T(1:length(Error_x_small_flux),4) = -10;
figure
subplot(1,2,1)
boxplot(T,'Labels',{'Big no flux','Big flux','Small no flux','Small flux'})
grid on
ylabel('[mm]')
ylim([0 4])
title('Error on x')
T = NaN(110,4);
% T(1:length(Error_y_big_noflux),1) = Error_y_big_noflux./MR_diam_big;
% T(1:length(Error_y_big_flux),2) = Error_y_big_flux./MR_diam_big;
% T(1:length(Error_y_small_noflux),3) = Error_y_small_noflux./MR_diam_small;
% T(1:length(Error_y_small_flux),4) = Error_y_small_flux./MR_diam_small;
T(1:length(Error_y_big_noflux),1) = Error_y_big_noflux;
T(1:length(Error_y_big_flux),2) = Error_y_big_flux;
T(1:length(Error_y_small_noflux),3) = Error_y_small_noflux;
T(1:length(Error_y_small_flux),4) = Error_y_small_flux;
subplot(1,2,2)
boxplot(T,'Labels',{'Big no flux','Big flux','Small no flux','Small flux'})
grid on
ylabel('[mm]')
ylim([0 4])
title('Error on y')

T = NaN(110,2);
T(1:length(Error_x_big_noflux),1) = Error_x_big_noflux;
T(1:length(Error_x_small_noflux),2) = Error_x_small_noflux;
figure
subplot(1,2,1)
boxplot(T,'Labels',{'Big no flux','Small no flux'})
grid on
ylabel('[mm]')
ylim([0 4])
title('Error on x')
T = NaN(110,2);
T(1:length(Error_y_big_noflux),1) = Error_y_big_noflux;
T(1:length(Error_y_small_noflux),2) = Error_y_small_noflux;
subplot(1,2,2)
boxplot(T,'Labels',{'Big no flux','Small no flux'})
grid on
ylabel('[mm]')
ylim([0 4])
title('Error on y')

%% Remove Outliers
% Error_x_big_noflux = rmoutliers(Error_x_big_noflux);
% Error_y_big_noflux = rmoutliers(Error_y_big_noflux);
% Error_x_big_flux = rmoutliers(Error_x_big_flux);
% Error_y_big_flux = rmoutliers(Error_y_big_flux);
% Error_x_small_noflux = rmoutliers(Error_x_small_noflux);
% Error_y_small_noflux = rmoutliers(Error_y_small_noflux);
% Error_x_small_flux = rmoutliers(Error_x_small_flux);
% Error_y_small_flux = rmoutliers(Error_y_small_flux);

%% ANOVA test
x_big_noflux = length(Error_x_big_noflux);
x_big_flux = length(Error_x_big_flux);
x_small_noflux = length(Error_x_small_noflux);
x_small_flux = length(Error_x_small_flux);
y_big_noflux = length(Error_y_big_noflux);
y_big_flux = length(Error_y_big_flux);
y_small_noflux = length(Error_y_small_noflux);
y_small_flux = length(Error_y_small_flux);
groups_x = [1*ones(1,x_big_noflux),2*ones(1,x_big_flux),3*ones(1,x_small_noflux),4*ones(1,x_small_flux)];
groups_y = [1*ones(1,y_big_noflux),2*ones(1,y_big_flux),3*ones(1,y_small_noflux),4*ones(1,y_small_flux)];
Error_x = [Error_x_big_noflux./MR_diam_big; Error_x_big_flux./MR_diam_big; Error_x_small_noflux./MR_diam_small; Error_x_small_flux./MR_diam_small];
Error_y = [Error_y_big_noflux./MR_diam_big, Error_y_big_flux./MR_diam_big, Error_y_small_noflux./MR_diam_small, Error_y_small_flux./MR_diam_small];
% Error_x = [Error_x_big_noflux; Error_x_big_flux; Error_x_small_noflux; Error_x_small_flux];
% Error_y = [Error_y_big_noflux, Error_y_big_flux, Error_y_small_noflux, Error_y_small_flux];
Error_big_noflux = [Error_x_big_noflux', Error_y_big_noflux];
groups_big_noflux = [1*ones(1,x_big_noflux),2*ones(1,y_big_noflux)];
Error_big_flux = [Error_x_big_flux',Error_y_big_flux];
groups_big_flux = [1*ones(1,x_big_flux),2*ones(1,y_big_flux)];
Error_small_noflux = [Error_x_small_noflux',Error_y_small_noflux];
groups_small_noflux = [1*ones(1,x_small_noflux),2*ones(1,y_small_noflux)];
Error_small_flux = [Error_x_small_flux',Error_y_small_flux];
groups_small_flux = [1*ones(1,x_small_flux),2*ones(1,y_small_flux)];

%ANOVA between x or y at different conditions
groups_x = groups_x>3;

p_x = anova1(Error_x,groups_x);
p_y = anova1(Error_y,groups_x);
%ANOVA between x and y fixing the conditions
% p_big_noflux = anova1(Error_big_noflux,groups_big_noflux);
% p_big_flux = anova1(Error_big_flux,groups_big_flux);
% p_small_noflux = anova1(Error_small_noflux,groups_small_noflux);
% p_small_flux = anova1(Error_small_flux,groups_small_flux);
groups_x = [1*ones(1,x_small_noflux),2*ones(1,x_small_flux)];
groups_y = [1*ones(1,y_small_noflux),2*ones(1,y_small_flux)];
Error_x = [Error_x_small_noflux; Error_x_small_flux];
Error_y = [Error_y_small_noflux, Error_y_small_flux];

Error = [Error_x';Error_y];
[D,P,stats] = manova1(Error',groups_x);
c1 = stats.canon(:,1);
c2 = stats.canon(:,2);
figure
gscatter(c2,c1,groups_x,'brg','o^s');
grid on
figure
manovacluster(stats)

Error2 = [groups_x',Error'];
maov1(Error2);


x_big_noflux = length(Error_x_big_noflux);
x_small_noflux = length(Error_x_small_noflux);
y_big_noflux = length(Error_y_big_noflux);
y_small_noflux = length(Error_y_small_noflux);
groups_x = [1*ones(1,x_big_noflux),2*ones(1,x_small_noflux)];
groups_y = [1*ones(1,y_big_noflux),2*ones(1,y_small_noflux)];
Error_x = [Error_x_big_noflux; Error_x_small_noflux];
Error_y = [Error_y_big_noflux, Error_y_small_noflux];
Error = [Error_x';Error_y];
[D2,P2,stats2] = manova1(Error',groups_x);


%% Statistics
data = Error_x_big_noflux;
mx = mean(data);
vx = var(data);
medx = median(data);
iqx = iqr(data);
disp(['The stats for BNF on x are: mean ',num2str(mx),'; var ',num2str(vx),'; median ',num2str(medx),'; IQ ',num2str(iqx)]);
data = Error_y_big_noflux;
my = mean(data);
vy = var(data);
medy = median(data);
iqy = iqr(data);
disp(['The stats for BNF on y are: mean ',num2str(my),'; var ',num2str(vy),'; median ',num2str(medy),'; IQ ',num2str(iqy)]);
data = Error_x_big_flux;
mx = mean(data);
vx = var(data);
medx = median(data);
iqx = iqr(data);
disp(['The stats for BF on x are: mean ',num2str(mx),'; var ',num2str(vx),'; median ',num2str(medx),'; IQ ',num2str(iqx)]);
data = Error_y_big_flux;
my = mean(data);
vy = var(data);
medy = median(data);
iqy = iqr(data);
disp(['The stats for BF on y are: mean ',num2str(my),'; var ',num2str(vy),'; median ',num2str(medy),'; IQ ',num2str(iqy)]);
data = Error_x_small_noflux;
mx = mean(data);
vx = var(data);
medx = median(data);
iqx = iqr(data);
disp(['The stats for SNF on x are: mean ',num2str(mx),'; var ',num2str(vx),'; median ',num2str(medx),'; IQ ',num2str(iqx)]);
data = Error_y_small_noflux;
my = mean(data);
vy = var(data);
medy = median(data);
iqy = iqr(data);
disp(['The stats for SNF on y are: mean ',num2str(my),'; var ',num2str(vy),'; median ',num2str(medy),'; IQ ',num2str(iqy)]);
data = Error_x_small_flux;
mx = mean(data);
vx = var(data);
medx = median(data);
iqx = iqr(data);
disp(['The stats for SF on x are: mean ',num2str(mx),'; var ',num2str(vx),'; median ',num2str(medx),'; IQ ',num2str(iqx)]);
data = Error_y_small_flux;
my = mean(data);
vy = var(data);
medy = median(data);
iqy = iqr(data);
disp(['The stats for SF on y are: mean ',num2str(my),'; var ',num2str(vy),'; median ',num2str(medy),'; IQ ',num2str(iqy)]);