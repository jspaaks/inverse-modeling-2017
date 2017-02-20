function visMatrixOfScatter(Sequences,ParRange,Extra,iGeneration)

nSequences = size(Sequences,3);
nPars = size(Sequences,2)-2;

nRowsSubplot = nPars-1;
nColsSubplot = nPars-1;

nGenerationsPlot = 250;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


defaultColors = [0,0,1;...
    0,0.5,0;
    1,0,0;...
    1,0.5,0;...
    1,0.9,0;...
    0.5882,0.7216,0.1804;...
    0,1,1;...
    1,0,1;...
    0.7529,0.7529,0.7529;...
    0.1490,0.5333,0.7176];

nDefaultColors = size(defaultColors,1);

if nSequences>nDefaultColors
    cVec = linspace(0.1,0.9,ceil((nSequences-nDefaultColors)^(1/3)));
    extraColors = allcomb(cVec,cVec,cVec);
else
    extraColors=[];
end

dreamSeqColors =[defaultColors;extraColors];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

set(gcf,'inverthardcopy','off',...
    'paperpositionmode','auto',...
    'numbertitle','off',...
    'name',['MatrixOfScatter (DREAM) - Figure ',num2str(gcf)])

canvasWidth = 0.9;
canvasHeight = 0.9;
canvasLeft = 0.1;
canvasBottom = 0.1;

axesWidth = canvasWidth/(nPars-1);
axesHeight = canvasHeight/(nPars-1);
axesWidthScaler = 0.95;
axesHeightScaler = 0.9;


clf



for iRowSubplot = 1:nRowsSubplot
    for iColSubplot = iRowSubplot:nColsSubplot

        L = canvasLeft+(iColSubplot-1)*axesWidth;
        B = canvasBottom+(nRowsSubplot-iRowSubplot)*axesHeight;
        W = axesWidth*axesWidthScaler;
        H = axesHeight*axesHeightScaler;

        axes('position',[L,B,W,H])
        
        for iSequence = 1:nSequences
            markerColor = dreamSeqColors(iSequence,1:3);
            plot(Sequences(max([1,iGeneration-nGenerationsPlot+1]):iGeneration,iColSubplot+1,iSequence),...
                Sequences(max([1,iGeneration-nGenerationsPlot+1]):iGeneration,iRowSubplot,iSequence),...
                'marker','s',...
                'markersize',4,...
                'linestyle','none',...
                'markerfacecolor',markerColor,...
                'markeredgecolor',markerColor);
            hold on
        end            
        

        set(gca,'xlim',[ParRange.minn(iColSubplot+1),ParRange.maxn(iColSubplot+1)],...
                'ylim',[ParRange.minn(iRowSubplot),ParRange.maxn(iRowSubplot)],...
                'color',[0.98,0.98,0.98])
            
        box on

        if iRowSubplot == iColSubplot
            if isfield(Extra,'ParNames')
                if isfield(Extra.ParNames,'interpreter')
                    xlabel(Extra.ParNames.labels{iColSubplot+1},...
                        'interpreter',Extra.ParNames.interpreter,...
                        'fontsize',14)
                    ylabel(Extra.ParNames.labels{iRowSubplot},...
                        'interpreter',Extra.ParNames.interpreter,...
                        'fontsize',14)
                else
                    xlabel(Extra.ParNames.labels{iColSubplot+1},...
                        'interpreter','none',...
                        'fontsize',14)
                    ylabel(Extra.ParNames.labels{iRowSubplot},...
                        'interpreter','none',...
                        'fontsize',14)
                end
            else
                xlabel(['\theta_{', num2str(iColSubplot+1),'}'],...
                    'interpreter','tex',...
                    'fontsize',14)
                ylabel(['\theta_{',num2str(iRowSubplot),'}'],...
                    'interpreter','tex',...
                    'fontsize',14)
            end
        else
            set(gca,'xticklabel',[],'yticklabel',[])
        end

            
    end
end

    