function [commModel]  = commonalityAnalysisNModels(megRDM,fmriRDM, targetRDM, nonTargetRDMs,validIndices)
%Apply commonality analysis 

%size(megRDM) = [durationMEG,widthRDM,widthRDM]
%size(fmriRDM) = [widthRDM,widthRDM]
%size(model1RDM) = [widthRDM,widthRDM]
%size(model2RDM) = [widthRDM,widthRDM]

%megRDM = 1
%fmriRDM = 2
%model1RDM = 3
%model2RDM = 4
%model3RDM = 5


durationMEG = size(megRDM,1);
widthRDM = size(megRDM,2);
commModel = zeros(durationMEG,1);


fmriRDM = fmriRDM(validIndices);
targetRDM = targetRDM(validIndices);
nonTargetRDMsFlat = zeros(size(nonTargetRDMs,1),length(fmriRDM));

for i = 1:size(nonTargetRDMs,1)
    thisNonTargetRDM = squeeze(nonTargetRDMs(i,:,:));
    nonTargetRDMsFlat(i,:) = thisNonTargetRDM(validIndices);
end
warning('off','all');
for i_time = 1:durationMEG
    thisMegRDM = squeeze(megRDM(i_time,:,:));
    
    
    Rsquared12N = multiRegressNvar(thisMegRDM(validIndices), ...
                                    [fmriRDM' ; nonTargetRDMsFlat]);
                                    
    Rsquared13N = multiRegressNvar(thisMegRDM(validIndices), ...
                                    [targetRDM' ; nonTargetRDMsFlat]);
                                                            
    if size(nonTargetRDMsFlat,1) >0
    Rsquared1N  = multiRegressNvar(thisMegRDM(validIndices), ...
                                    nonTargetRDMsFlat);
    else
        Rsquared1N = 0;
    end
                                
    Rsquared123N = multiRegressNvar(thisMegRDM(validIndices), ...
                                    [fmriRDM' ; targetRDM' ; nonTargetRDMsFlat]);
                                    
    
    commModel(i_time) = Rsquared12N + Rsquared13N  - Rsquared1N - Rsquared123N;
        
end
warning('on','all');