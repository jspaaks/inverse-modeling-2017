
restoredefaultpath
addpath('.\..\dream')

clear
close all
clc

load('leafriver.mat')

iStart = 66
iEnd = iStart+181


MCMCPar.n = 5;                          % Dimension of the problem
MCMCPar.seq = 10;                       % Number of Markov Chains / sequences
MCMCPar.ndraw = MCMCPar.n*600;          % Maximum number of function evaluations
MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
MCMCPar.DEpairs = 3;                    % Number of DEpairs
MCMCPar.steps = 10;                     % Number of steps in sem
MCMCPar.eps = 5e-2;                     % Random error for ergodicity
MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

% Give the parameter ranges (minimum and maximum values)
ParRange.minn = [  1.0, 0.1, 0.10, 0.000, 0.100];
ParRange.maxn = [500.0, 2.0, 0.99, 0.100, 0.990];

% Define the measured data
Measurement.MeasData = dailyDischarge(iStart:iEnd,1); 
Measurement.N = size(Measurement.MeasData,1);

% Define modelName
ModelName = 'hymod';
Extra.dailyDischarge = dailyDischarge;
Extra.dailyPotEvapTrans = dailyPotEvapTrans;
Extra.dailyPrecip = dailyPrecip;
Extra.iStart = iStart;
Extra.iEnd = iEnd;
Extra.numTime = numTime;

Extra.BoundHandling = 'Reflect';        % Define the boundary handling
Extra.save_in_memory = 'Yes';
Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
Extra.InitPopulation = 'LHS_BASED';     % What type of initial sampling


% Define likelihood function -- Sum of Squared Error
option = 3; 
Measurement.Sigma = [];

% specify the visualization routine
Extra.visScriptName = 'visDream';
% Specify labels for the parameter names
Extra.ParNames.labels = {'c_{max}','b_{exp}','\alpha_{q}','R_s','R_q'};
Extra.ParNames.interpreter = 'tex';

% Scale of Delayed Rejection if used in algorithm
Extra.DR = 'No'; 
Extra.DRscale = 10; 
Restart = 'No';

% Run the MCMC algorithm
[Sequences,Reduced_Seq,X,output,hist_logp] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option,Restart);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % %                 ANALYSIS AND VISUALIZATION                  % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

parCols=1:MCMCPar.n;
lastRows = [-99:0]+size(Sequences,1);
parSets = [];
for iSeq=1:MCMCPar.seq
    parSets = cat(1,parSets,Sequences(lastRows,parCols,iSeq));
end

for k=1:size(parSets,1)
    parVec = parSets(k,:);
    ySim(k,:) = hymod(parVec,Extra);
end


Y = prctile(ySim,[2.5,50,97.5],1);

figure(123)
subplot(2,1,1)
plot(numTime(iStart:iEnd),dailyPotEvapTrans(iStart:iEnd),'-r',...
     numTime(iStart:iEnd),0.1*dailyPrecip(iStart:iEnd),'-c' )
hLeg = legend('dailyPotEvapTrans [mm\cdot{}day^{-1}]','dailyPrecip [10\cdot{}mm\cdotday^{-1}]');
set(hLeg,'fontsize',8,'color','none','xcolor','w','ycolor','w','location','northwest')

title('Leaf River watershed observations')    
datetick('x')

subplot(2,1,2)
h=plot(numTime(iStart:iEnd),dailyDischarge(iStart:iEnd),'om',...
       numTime(iStart:iEnd),Y(1,:),'--k',...
       numTime(iStart:iEnd),Y(2,:),'-k',...
       numTime(iStart:iEnd),Y(3,:),'--k');

set(gca,'ylim',[0,420])   
set(h(1),'markersize',3,'markerfacecolor','m')
ylabel('dailyDischarge [m^3\cdot{}s^{-1}]')
hLeg = legend(h([1,2,3]),'observed discharge','95% parameter uncertainty','median prediction');
set(hLeg,'fontsize',8,'color','none','xcolor','w','ycolor','w','location','northwest')
xlabel('time')

datetick('x')
drawnow




figure(123)
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print('prediction.png','-dpng','-r300')

figure(2)
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print('plotSeq.png','-dpng','-r300')

figure(3)
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print('matrixOfScatterParPar.png','-dpng','-r300')
