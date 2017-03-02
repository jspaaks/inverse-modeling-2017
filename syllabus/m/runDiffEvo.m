% differential evolution script

% restore the default MATLAB path in case it was changed:
restoredefaultpath

% add the folder to the MATLAB path that contains the 
% data needed for some of the benchmark functions:
addpath(['.',filesep,'benchmark-data'])

% clear any old variables:
clear

% close any old figures:
close all

% clen up the command window, so we don't get 
% confused by any old error messages:
clc


% the benchmark function needs a few variables:
global initial_flag; 
initial_flag = 0;           

% this next command sets the flag, indicating which function is used:
funcFlag = 6;

nGenerations = 2;
nPop = 10;
nDims = 2;
parents = repmat(NaN, [nPop, nDims + 1]);
parCols = 1:nDims;
objCol = nDims + 1;

parSpaceLowerBounds = [-100, -100];
parSpaceUpperBounds = [100, 100];
parSpaceRange = parSpaceUpperBounds - parSpaceLowerBounds;
F = 0.6;
K = 0.4;

parents(1:nPop,1:nDims) = repmat(parSpaceLowerBounds,[nPop, 1]) + rand(nPop,nDims) .* repmat(parSpaceRange,[nPop, 1]);


for iPop = 1:nPop
    parVec = parents(iPop,1:nDims);
    objScore = benchmark_func(parVec,funcFlag);
    parents(iPop, objCol) = objScore;
end

for iGeneration = 1:nGenerations
    
    proposals = repmat(NaN, [nPop, nDims + 1]);

    for iPop = 1:nPop

        v = randperm(nPop, 3);
        r1 = parents(v(1), parCols);
        r2 = parents(v(2), parCols);
        r3 = parents(v(3), parCols);

        dist1 = r1 - parents(iPop, parCols);
        dist2 = r3 - r2;

        proposals(iPop, parCols) = parents(iPop, parCols) + F * dist1 + K * dist2;
    end

    for iPop = 1:nPop
        parVec = proposals(iPop,1:nDims);
        objScore = benchmark_func(parVec,funcFlag);
        proposals(iPop, objCol) = objScore;
    end

    children = repmat(NaN, [nPop, nDims + 1]); 
    for iPop = 1:nPop
        if proposals(iPop,objCol) < parents(iPop,objCol)
            children(iPop, parCols) = proposals(iPop, parCols);
            children(iPop, objCol) = proposals(iPop, objCol);
        else
            children(iPop, parCols) = parents(iPop, parCols);
            children(iPop, objCol) = parents(iPop, objCol);
        end
    end

    if iGeneration == 1
        figure(111)
        set(gcf, 'Position', [681, 304, 1000, 300])
        clf
        nSubplots = 4;
        for iSubplot = 1:nSubplots-1
            W = 1 / (nSubplots + 1) - 0.075;
            L = 0.075 + (iSubplot - 1) * 1 / (nSubplots + 1);
            B = 0.15;
            H = 0.70;
            axes('position', [L,B,W,H]);
            hold on
            if ismember(iSubplot,[1,2,3])
                plot(parents(:,1),parents(:,2), '+k', 'MarkerSize',6);
            end
            if ismember(iSubplot,[2,3])
                plot(proposals(:,1),proposals(:,2), 'sk', 'MarkerSize',4, 'MarkerFaceColor', 'k');
                plot([parents(:,1),proposals(:,1)]',[parents(:,2),proposals(:,2)]', 'LineStyle','-','Color',[0.5,0.5,0.5]);
            end
            if ismember(iSubplot,[3])
                plot(children(:,1), children(:,2),'ok', 'MarkerSize',9);
            end
            if iSubplot == 1
                text(-160,110,'\bf{A}')
            elseif iSubplot == 2
                text(-160,110,'\bf{B}')
            elseif iSubplot == 3
                text(-160,110,'\bf{C}')
            end
            box on
            xlabel('par_{1}')
            ylabel('par_{2}')
            title(['iGeneration = ', num2str(iGeneration)])
            set(gca, 'XLim', [parSpaceLowerBounds(1), parSpaceUpperBounds(1)])
            set(gca, 'YLim', [parSpaceLowerBounds(2), parSpaceUpperBounds(2)])
        end
    else
        iSubplot = iSubplot + 1;
        L = 0.075 + (iSubplot - 1) * 1 / (nSubplots + 1);
        axes('position', [L,B,W,H]);
        plot(parents(:,1),parents(:,2), '+k', 'MarkerSize',6);
        xlabel('par_{1}')
        ylabel('par_{2}')
        title(['iGeneration = ', num2str(iGeneration)])
        set(gca, 'XLim', [parSpaceLowerBounds(1), parSpaceUpperBounds(1)])
        set(gca, 'YLim', [parSpaceLowerBounds(2), parSpaceUpperBounds(2)])
        text(-160,110,'\bf{D}')
        box on

        iSubplot = iSubplot + 1;
        L = 0.075 + (iSubplot - 1) * 1 / (nSubplots + 1) - 0.05;
        axes('position', [L,B,W,H]);
        plot(0.5, 1, '+k', 'MarkerSize',6)
        hold on
        text(1.0, 1, 'parents')
        plot(0.5, 2, 'sk', 'MarkerSize',4, 'MarkerFaceColor', 'k');
        text(1.0, 2, 'proposals')
        plot([0.25,0.75], [3,3], 'LineStyle','-','Color',[0.5,0.5,0.5]);
        text(1.0, 3, ['parents and their', char(10), 'associated ',char(10),'proposals'])
        plot(0.5,4,'ok', 'MarkerSize',9);
        text(1.0, 4, 'accepted children')
        set(gca, 'XLim', [0,5], 'YLim', [0,5], 'YDir', 'reverse')
        axis off
        set(gcf,'paperpositionmode','auto','inverthardcopy','off')
        print('diff-evo-principle.eps','-depsc2','-r300','-loose')
        
    end

    parents = children;
end


