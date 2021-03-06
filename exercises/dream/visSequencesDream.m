function visSequencesDream(Sequences,ParRange,Extra,iGeneration)

nSequences = size(Sequences,3);
nPars = size(Sequences,2)-2;

nRowsSubplot = min([4,nPars]);
nColsSubplot = ceil(nPars/nRowsSubplot);


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
    'name',['Sequences (DREAM) - Figure ',num2str(gcf)])

clf

for iColSubplot = 1:nColsSubplot
    for iRowSubplot = 1:nRowsSubplot

        iPar = (iColSubplot-1)*nRowsSubplot+iRowSubplot;
        if iPar>nPars
            continue
        end
            
        iSubplot = (iRowSubplot-1)*nColsSubplot+iColSubplot;

        subplot(nRowsSubplot,nColsSubplot,iSubplot)
        for iSequence = 1:nSequences
            markerColor = dreamSeqColors(iSequence,1:3);
            xSeq = iSequence+0:nSequences:iGeneration*nSequences;
            plot(xSeq,Sequences(1:iGeneration,iPar,iSequence),...
                'marker','o',...
                'markersize',3,...
                'linestyle','none',...
                'markerfacecolor',markerColor,...
                'markeredgecolor',markerColor);
            hold on
        end

        set(gca,'ylim',[ParRange.minn(iPar),ParRange.maxn(iPar)],...
            'color',[0.98,0.98,0.98])

        if isfield(Extra,'ParNames')
            if isfield(Extra.ParNames,'interpreter')
                ylabel(Extra.ParNames.labels(iPar),'interpreter',Extra.ParNames.interpreter)
            else
                ylabel(Extra.ParNames.labels(iPar),'interpreter','none')
            end
        else
            ylabel(['\theta_{',num2str(iPar),'}'],'interpreter','tex')
        end


        if iRowSubplot<nRowsSubplot & iPar<nPars
            set(gca,'xticklabel','')
        else
            xlabel('number of function evaluations')
        end
    end
end

