function visObjScoreDream(Sequences,iGeneration)

clf
nSequences = size(Sequences,3);
nPars = size(Sequences,2)-2;


% if isequal(log(Sequences(1:iGeneration,nPars+1,:)),Sequences(1:iGeneration,nPars+2,:))
%     subplotVec = [1,2];
% elseif isequal(Sequences(1:iGeneration,nPars+1,:),Sequences(1:iGeneration,nPars+2,:))
    subplotVec = 1;
% else
%     error('uhoh')
% end

for iSubplot = subplotVec
    
    subplot(subplotVec(end),1,iSubplot)
    objCol = nPars+iSubplot;

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
        'name',['Objective scores (DREAM) - Figure ',num2str(gcf)])


    markerColor = [0.1,0.1,1];
    xObjScore = [1:iGeneration*nSequences]';
    yObjScoreTMP = shiftdim(Sequences(1:iGeneration,objCol,:),2);

    for iSequence = 1:nSequences
        markerColor = dreamSeqColors(iSequence,1:3);
        xSeq = iSequence+0:nSequences:iGeneration*nSequences;
        ySeq = Sequences(1:iGeneration,objCol,iSequence)';
        plot(xSeq,ySeq,...
            'marker','o',...
            'markersize',3,...
            'linestyle','none',...
            'markerfacecolor',markerColor,...
            'markeredgecolor',markerColor);
        hold on
    end

    set(gca,'color',[0.98,0.98,0.98],'yscale','log')
    %set(gca,'color',[0.98,0.98,0.98],'yscale','linear')
    xlabel('number of function evaluations')
    
    if iSubplot==1
        ylabel('objScore')
    elseif iSubplot==2
        ylabel('log(objScore)')
    else
    end
    
end

        







