

% user-settable variables:
visualizeObjScores = true;
visualizeSequences = true;
visualizeMatrixOfScatter = true;
visualizeStateSpace = true;



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


if exist('establishedFigures','var')==0
    
    if visualizeObjScores
        scemuaSubplotScreen(2,3,1,'figureNumber',1)
    end
    if visualizeStateSpace
        scemuaSubplotScreen(2,3,4,'figureNumber',123)
    end
    
    if visualizeSequences
        scemuaSubplotScreen(1,3,2,'figureNumber',2)
    end
    if visualizeMatrixOfScatter
        scemuaSubplotScreen(1,3,3,'figureNumber',3)
    end
    
    establishedFigures = true;
end
    

%iGeneration = counter;
iGeneration = find(all(all(Sequences==0,2),3),1,'first')-1;
if isempty(iGeneration)
    iGeneration = size(Sequences,1);
end

% if counter~=iGeneration
%     error('something''s worng here')
% end

if exist('iGenerationLastVis','var')==0
    iGenerationLastVis = 0;
end
if exist('visInterval','var')==0
    visInterval = 50;  % measured in generations, not function evaluations
end

if (iGeneration-iGenerationLastVis)>=visInterval || iGeneration==size(Sequences,1)
    if strcmp(Extra.save_in_memory,'Yes')
        
        if visualizeObjScores
            figure(1)
            visObjScoreDream(Sequences,iGeneration)
        end


        if visualizeSequences
            figure(2)
            visSequencesDream(Sequences,ParRange,Extra,iGeneration)
        end
        
        if visualizeMatrixOfScatter
            figure(3)
            visMatrixOfScatter(Sequences,ParRange,Extra,iGeneration)
        end
        
        iGenerationLastVis = iGeneration;
 
        drawnow
    else
        error('No information to visualize')
    end
end        

if iGeneration==size(Sequences,1)
    clear establishedFigures
end
