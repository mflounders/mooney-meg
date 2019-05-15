function [commModel, fMRIMEG]  = combinedNModels(megRDM,fmriRDM, modelsRDM,validIndices)
%Calculates combined variance for given RDMs

%size(megRDM) = [durationMEG,widthRDM,widthRDM]
%size(fmriRDM) = [widthRDM,widthRDM]
%size(model1RDM) = [widthRDM,widthRDM]
%size(model2RDM) = [widthRDM,widthRDM]

durationMEG = size(megRDM,1);
widthRDM = size(megRDM,2);
lowerTriangleIndices = validIndices; %indices for lower triangle of RDM

fmriRDM = fmriRDM(lowerTriangleIndices);
modelsRDMsFlat = zeros(size(modelsRDM,1),length(fmriRDM));
for i = 1:size(modelsRDM,1)
    thisModelRDM = squeeze(modelsRDM(i,:,:));
    modelsRDMsFlat(i,:) = thisModelRDM(lowerTriangleIndices);
end

commModel = zeros(durationMEG,1);
fMRIMEG = zeros(durationMEG,1);

for i_time = 1:durationMEG
    %i_time
    thisMegRDM = squeeze(megRDM(i_time,:,:));
     
    Rsquared12 = multiRegressNvar(thisMegRDM(lowerTriangleIndices), ...
                                    fmriRDM');
    
    Rsquared1N = multiRegressNvar(thisMegRDM(lowerTriangleIndices), ...
                                    modelsRDMsFlat);
    
    Rsquared12N = multiRegressNvar(thisMegRDM(lowerTriangleIndices), ...
                                    [fmriRDM' ; modelsRDMsFlat]);
    
    commModel(i_time) = Rsquared12 + Rsquared1N - Rsquared12N;
    fMRIMEG(i_time) = Rsquared12;
end
