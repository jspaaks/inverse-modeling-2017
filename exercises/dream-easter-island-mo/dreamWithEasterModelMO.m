
restoredefaultpath
addpath('.\..\dream')

clear
close all
clc

dataFileRes = 'pop-res-data-ref.mat';
dataFilePop = 'pop-res-data-high-deathrate.mat';

load(dataFilePop,'PopMeas')
load(dataFileRes,'ResourceMeas','tMeas')



MCMCPar.n = 6;                          % Dimension of the problem
MCMCPar.seq = 10;                       % Number of Markov Chains / sequences
MCMCPar.ndraw = MCMCPar.n*300;         % Maximum number of function evaluations
MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
MCMCPar.DEpairs = 3;                    % Number of DEpairs
MCMCPar.steps = 10;                     % Number of steps in sem
MCMCPar.eps = 5e-2;                     % Random error for ergodicity
MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes


% Give the parameter ranges (minimum and maximum values)
ParRange.minn = [0.02,0.03,160,150,0.004,1e-7]/3;
ParRange.maxn = [0.02,0.03,160,150,0.004,1e-7]*3;

% Define the measured data
Measurement.MeasData = PopMeas; 
Measurement.N = size(Measurement.MeasData,1);
Measurement.Sigma = [];

% Define modelName
ModelName = 'eastermodelMO';
Extra.tMeas = tMeas;
Extra.tStart = 0;
Extra.tEnd = 2000;
Extra.dt = 1.0;
Extra.PopMeas = PopMeas;
Extra.ResMeas = ResourceMeas;
Extra.dontPlot = false;

Extra.BoundHandling = 'Reflect';        % Define the boundary handling
Extra.save_in_memory = 'Yes';
Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
Extra.InitPopulation = 'LHS_BASED';     % What type of initial sampling

% Define likelihood function -- this is now a direct probability,
% because I'm calculating the SSE inside eastermodelMO.
option = 3; 

% specify the visualization routine
Extra.visScriptName = 'visDream';
% Specify labels for the parameter names
Extra.ParNames.labels = {'rel. birth rate','rel. death rate',...
    'carrying capacity','consumption response','rel. growth rate',...
    'fitting parameter'};
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
if ~Extra.argOutIsErr

    parCols=1:MCMCPar.n;
    lastRows = [-99:0]+size(Sequences,1);
    parSets = [];
    for iSeq=1:MCMCPar.seq
        parSets = cat(1,parSets,Sequences(lastRows,parCols,iSeq));
    end

    Extra.dontPlot = true;
    u=1;
    for k=1:size(parSets,1)
        parVec = parSets(k,:);
        TMP = eastermodelMO(parVec,Extra);
        if all(TMP==-999)
            %
        else
            ySimPop(u,:) = TMP(:,1);
            ySimRes(u,:) = TMP(:,2);
            u = u + 1;
        end
    end


    probPop = prctile(ySimPop,[2.5,50,97.5],1);
    probRes = prctile(ySimRes,[2.5,50,97.5],1);

    tSim = Extra.tStart:Extra.dt:Extra.tEnd;

    figure(123)
    clf
    subplot(2,1,1)
    measStyle = {'marker','o','markersize',3,'linestyle','none'};
    hResObs = plot(tMeas,ResourceMeas,'markeredgecolor',[0,0.5,0],...
        'markerfacecolor',[0,0.5,0],measStyle{:});
    hold on
    h2 = plot(tSim(1,1:end-1),probRes(1,:),'--k',...
              tSim(1,1:end-1),probRes(2,:),'-k',...
              tSim(1,1:end-1),probRes(3,:),'--k');
    set(gca,'xlim',[0,2000],'ylim',[0,250])
    hLegR = legend([hResObs;h2(1:2)],'ResourceMeas','95% parameter uncertainty','median prediction');
    set(hLegR,'color','none','location','NorthWest','xcolor','w','ycolor','w','fontsize',8);
    xlabel('time [years]')
    ylabel('Biomass')
    title(dataFileRes)


    subplot(2,1,2)
    hPopObs = plot(tMeas,PopMeas,'markeredgecolor',[1,0,1],...
        'markerfacecolor',[1,0,1],measStyle{:});
    hold on
    h3 = plot(tSim(1,1:end-1),probPop(1,:),'--k',...
              tSim(1,1:end-1),probPop(2,:),'-k',...
              tSim(1,1:end-1),probPop(3,:),'--k');
    set(gca,'xlim',[0,2000],'ylim',[0,15000])
    hLegP = legend([hPopObs;h3(1:2)],'PopMeas','95% parameter uncertainty','median prediction');
    set(hLegP,'color','none','location','NorthWest','xcolor','w','ycolor','w','fontsize',8);
    xlabel('time [years]')
    ylabel('Population')
    title(dataFilePop)



    figure(1)
    set(gcf,'paperpositionmode','auto','inverthardcopy','off')
    print(['objscores.png'],'-dpng','-r300')

    figure(2)
    set(gcf,'paperpositionmode','auto','inverthardcopy','off')
    print(['sequences.png'],'-dpng','-r300')

    figure(3)
    set(gcf,'paperpositionmode','auto','inverthardcopy','off')
    print(['matrixOfScatter.png'],'-dpng','-r300')


    figure(123)
    set(gcf,'paperpositionmode','auto','inverthardcopy','off')
    print('pop-res-figures.png','-dpng','-r300')
    print('pop-res-figures.eps','-depsc2','-r300','-loose')

end


