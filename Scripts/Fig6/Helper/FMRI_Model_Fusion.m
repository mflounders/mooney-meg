function FMRI_Model_Fusion
%Correlates fMRI RDMs with the three models.

load('Data/Raw/MEG_RDM.mat')
load('Data/Raw/FMRI_RDM.mat')
load(['Data/Raw/Models.mat'])
lowerTriangleInds = find(tril(ones(99),-1));
numROIs = 20;
numModels = 3;
rFMRIModel = zeros(numModels,numROIs);
pFMRIModel = zeros(numModels,numROIs);
for i = 1:3
    thisModel = squeeze(models(i,:,:));
    validInds = lowerTriangleInds(find(~isnan(thisModel(lowerTriangleInds))));
    for j = 1:20
        thisFmriROI = squeeze(FMRI_RDM(:,:,j));
        [rFMRIModel(i,j), pFMRIModel(i,j)] = corr(thisModel(validInds),thisFmriROI(validInds),'type','spearman');
    end
end
        
        
save(['Data/Processed/Figure6/FRMI_Model_Fusion.mat'],'rFMRIModel','pFMRIModel');