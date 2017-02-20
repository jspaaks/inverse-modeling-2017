
restoredefaultpath
addpath('.\..\dream')

clear
close all
clc

% create the artificial data:
a = 0.43
b = 9.87

xObs = [-3:0.1:7]'
yObsTrue = a*xObs+b
gaussTerm = randn(size(yObsTrue))*0.3
yObs = yObsTrue + gaussTerm


% visualize the artificial data:
figure(123)
h1=plot(xObs,yObsTrue,'-b',xObs,yObs,'.m');
xlabel('x')
ylabel('y')

% DREAM algorithm settings:
MCMCPar.n = 2;                          % Dimension of the problem
MCMCPar.seq = 10;                       % Number of Markov Chains / sequences
MCMCPar.ndraw = MCMCPar.n*2000;         % Maximum number of function evaluations
MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
MCMCPar.DEpairs = 3;                    % Number of DEpairs
MCMCPar.steps = 10;                     % Number of steps in sem
MCMCPar.eps = 5e-2;                     % Random error for ergodicity
MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes


% Give the parameter ranges (minimum and maximum values)
ParRange.minn = [-3.0,-10.0];
ParRange.maxn = [ 3.0, 20.0];

% Define the observations:
Measurement.MeasData = yObs; 
Measurement.N = size(Measurement.MeasData,1);
Measurement.Sigma = [];

% Define likelihood function -- Sum of Squared Error
option = 3; 

% Define modelName
ModelName = 'linearmodel';
Extra.xObs = xObs;



Extra.BoundHandling = 'Reflect';        % Define the boundary handling
Extra.save_in_memory = 'Yes';
Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
Extra.InitPopulation = 'LHS_BASED';     % What type of initial sampling


% specify the visualization routine
Extra.visScriptName = 'visDream';
% Specify labels for the parameter names
Extra.ParNames.labels = {'slope','intercept'};
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
    ySim(k,:) = linearmodel(parVec,Extra);
end


Y = prctile(ySim,[2.5,50,97.5],1);
figure(123)
hold on
h2 = plot(xObs,Y(1,:),'--k',...
     xObs,Y(2,:),'-k',...
     xObs,Y(3,:),'--k');
hLeg = legend([h1;h2(1:2)],'true','observed','95% parameter uncertainty','median prediction');
set(hLeg,'color','none','box','off','location','NorthWest');





