
restoredefaultpath
addpath('.\..\dream')

clear
close all
clc

% load('bound_ideal.mat')
% bound = bound_ideal;
% clear bound_ideal
% load('stor_ideal.mat')
% stor = stor_ideal;
% clear stor_ideal

% load('bound_211_short.mat')
% bound = bound_211;
% clear bound_211
% load('stor_211_short.mat')
% stor = stor_211;
% clear stor_211

load('bound_211.mat')
bound = bound_211;
clear bound_211
load('stor_211.mat')
stor = stor_211;
clear stor_211

% load('bound_210.mat')
% bound = bound_210;
% clear bound_210
% load('stor_210.mat')
% stor = stor_210;
% clear stor_210


meanTime = mean([stor(1:end-1,1),stor(2:end,1)],2);
storDeriv = (stor(2:end,2)-stor(1:end-1,2))./(stor(2:end,1)-stor(1:end-1,1));


yMax = 3.5;
scaledPrec = (bound(:,2)-min(bound(:,2)))./...
    (max(bound(:,2))-min(bound(:,2)));
scaledEvap = (bound(:,3)-min(bound(:,3)))./...
    (max(bound(:,3))-min(bound(:,3)));

figure(123)
h1 = plot(bound(:,1),yMax-0.9*scaledPrec*yMax,'-c.',...
          bound(:,1),0.9*scaledEvap*yMax,'-r.',...
          stor(:,1),stor(:,2),'mo');
set(h1(3),'markersize',4,'markerfacecolor','w')
title('interception model')



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
ParRange.minn = [0.4,   0.0,0.0,0.0];
ParRange.maxn = [1.0,1000.0,7.0,1.0];

% Define the measured data
Measurement.MeasData = storDeriv(:,1); 
Measurement.N = size(Measurement.MeasData,1);

% Define modelName
ModelName = 'interceptionmodel';
Extra.bound = bound;
Extra.stor = stor;
Extra.argOutIsDeriv = true;


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
Extra.ParNames.labels = {'vegFraction','drainEff','canopyMax','evapFactor'};
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


Extra.argOutIsDeriv = false;
for k=1:size(parSets,1)
    parVec = parSets(k,:);
    ySim(k,:) = interceptionmodel(parVec,Extra);
end


Y = prctile(ySim,[2.5,50,97.5],1);
figure(123)
hold on
h2 = plot(bound(:,1),Y(1,:),'--k',...
          bound(:,1),Y(2,:),'-k',...
          bound(:,1),Y(3,:),'--k');
hLeg = legend([h1;h2(1:2)],'precipitation (scaled) ','pot. evap. (scaled)',...
    'observed storage','95% parameter uncertainty','median prediction');
set(hLeg,'location','NorthWest','fontsize',8,'color','none','xcolor','w','ycolor','w');
xlabel('time [days]')
ylabel('canopy storage [mm]')





figure(123)
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print('prediction.png','-dpng','-r300')

figure(2)
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print('plotSeq.png','-dpng','-r300')

figure(3)
set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print('matrixOfScatterParPar.png','-dpng','-r300')
