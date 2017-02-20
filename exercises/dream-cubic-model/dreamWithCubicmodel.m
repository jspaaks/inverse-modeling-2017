
restoredefaultpath
addpath('.\..\dream')

clear
close all
clc

% true parameters:
a = -1;
b = 7;
c = -5;
d = 9;

% you only need to change xObsStart, xStep, xEnd, and randFactor:
randnFactor = 5;
xObsStart = -1;
xObsStep = 0.6;
xObsEnd = 7;



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % %           no need to change anything below this line            % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


xTrue = [min([xObsStart,-1]):0.01:max([xObsEnd,7])]';
yTrue = a*xTrue.^3+b*xTrue.^2+c*xTrue+d;
randn('state',0)
randnLarge = randn(size(xTrue));

xObs = [xObsStart:xObsStep:xObsEnd]';
yObs = a*xObs.^3+b*xObs.^2+c*xObs+d;

% round off to 0.001 otherwise the ismember trick fails sometimes:
xTrue = round(xTrue*1000)/1000;
xObs = round(xObs*1000)/1000;

% make sure the rand values are consistent between exercises:
gaussTerm = randnLarge(ismember(xTrue,xObs))*randnFactor;

yObsNoisy = yObs + gaussTerm;


figure(123)
clf
h1 = plot(xTrue,yTrue,'-b',...
          xObs,yObs,'.m',...
          xObs,yObsNoisy,'om');

set(gca,'xlim',[xTrue(1),xTrue(end)],'ylim',[-20,70])

title('cubic function: y=ax^3+bx^2+cx+d')


MCMCPar.n = 4;                          % Dimension of the problem
MCMCPar.seq = 10;                       % Number of Markov Chains / sequences
MCMCPar.ndraw = MCMCPar.n*1000;         % Maximum number of function evaluations
MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
MCMCPar.DEpairs = 3;                    % Number of DEpairs
MCMCPar.steps = 10;                     % Number of steps in sem
MCMCPar.eps = 5e-2;                     % Random error for ergodicity
MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes


% Give the parameter ranges (minimum and maximum values)
ParRange.minn = [-20.0,-40.0,-80.0,-120.0];
ParRange.maxn = [ 20.0, 40.0, 80.0, 120.0];

% Define the measured streamflow data
Measurement.MeasData = yObsNoisy; 
Measurement.N = size(Measurement.MeasData,1);

% Define modelName
ModelName = 'cubicmodel';
Extra.xObs = xObs;

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
Extra.ParNames.labels = {'a','b','c','d'};
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

Extra.xObs = xTrue;
clear ySim
for k=1:size(parSets,1)
    parVec = parSets(k,:);
    ySim(k,:) = cubicmodel(parVec,Extra);
end


Y = prctile(ySim,[2.5,50,97.5],1);
figure(123)
hold on
h2 = plot(Extra.xObs,Y(1,:),'--k',...
          Extra.xObs,Y(2,:),'-k',...
          Extra.xObs,Y(3,:),'--k');
xlabel('x')
ylabel('y')
set(gca,'xlim',[xTrue(1),xTrue(end)],'yLim',[-20,70]);
hLeg = legend([h1;h2(1:2)],'true','observed, no noise',...
    'observed, noisy','95% parameter uncertainty',...
'median prediction');
set(hLeg,'color','none','box','off','location','NorthWest','fontsize',8);

figure(123)
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print(['model-prediction.png'],'-dpng','-r300')
figure(3)
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print(['matrix-of-scatter.png'],'-dpng','-r300')




