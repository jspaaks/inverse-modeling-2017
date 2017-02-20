clear
clc
close all

warning off;

% set the random seed
rand('state',0)

% Load the Leaf River data
load bound.txt;

% Then read the boundary conditions
Extra.PET = bound(:,5);
Extra.Precip = sum(bound(:,6:9),2);

% Define the measured streamflow data
Extra.MeasData = bound(:,4);

Extra.numTime = datenum(bound(:,3),bound(:,2),bound(:,1));
strTime = datestr(Extra.numTime);

figure
subplot(3,1,1)
plot(Extra.numTime,bound(:,5),'-r')
ylabel('PET')

subplot(3,1,2)
plot(Extra.numTime,sum(bound(:,6:9),2),'-c')
ylabel('Precip')

subplot(3,1,3)
plot(Extra.numTime,bound(:,4),'-b')
set(gca,'yscale','log')
ylabel('Q')

% Give the parameter ranges (minimum and maximum values)
ParRange.minn = [1.0 0.10 0.10 0.00 0.10];
ParRange.maxn = [500 2.00 0.99 0.10 0.99];

% Define the subsets
Ix_jan01 = find(bound(:,1)==1&bound(:,2)==1&bound(:,3)==1953);
Ix_mar31 = find(bound(:,1)==31&bound(:,2)==3&bound(:,3)==1953);

subset = [Ix_jan01,Ix_mar31]

% Select the first subset for the multistart Simplex
startRow = subset(1,1);
endRow = subset(1,2);

disp(['Using data from ',strTime(startRow,:),' to ',strTime(endRow,:),'.'])

Extra.MaxT = endRow;
Extra.calPeriod = [startRow,endRow];

% Select 20 start points and run Simplex

for i=1:20
    ini_par=rand([1 5]).*(ParRange.maxn-ParRange.minn)+ParRange.minn;
    options = optimset('display','final', 'MaxIter',5000, 'MaxFunEvals', 15000);
    [X(i,1:5),FVAL(i)]= fminsearch('hymod_OF',ini_par,options,Extra);
end

X
FVAL

% plot the 20 simulations
figure
for i=1:20
    [SimData(i,:)]=hymod(X(i,:),Extra);
end
plot(Extra.calPeriod(1):Extra.calPeriod(2),Extra.MeasData(Extra.calPeriod(1):Extra.calPeriod(2)),'.m',...
    Extra.calPeriod(1):Extra.calPeriod(2),SimData);
ylabel('Discharge [m^3/s]')
xlabel('Calibration Period [days]')
title('dots=measured, lines=simulations')
% plot the 20 parameter sets

figure
minn = repmat(ParRange.minn,[size(X,1),1]);
maxn = repmat(ParRange.maxn,[size(X,1),1]);
plot(1:5,(X-minn)./(maxn-minn),'.-')
set(gca,'xtick',[1:5],'xticklabel',{'cmax','bexp','alpha','Rs','Rq'})
ylabel('Optimized model parameters')
xlabel('parameter names')


