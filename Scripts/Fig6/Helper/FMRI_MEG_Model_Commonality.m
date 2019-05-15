function FMRI_MEG_Model_Commonality(i_ROI)

%Applies Commmonality analysis to FMRI, MEG, and Model RDMs

%% PRESETS
lengthTrial = 550;
picked3Models = [1 2 3];
numPicked3Models = length(picked3Models);
picked2Models = [1 3];
numPicked2Models = length(picked2Models);
numPermutations = 1000;
widthRDM = 99;

%% Load Data
load('Data/Raw/MEG_RDM.mat');
load('Data/Raw/MEG_Permutations.mat');   %Contains shuffled indices for each permutation
load('Data/Raw/FMRI_RDM.mat');
load(['Data/Raw/Models.mat']);

groupAverage_MEG_RDM = mean(MEG_RDM,4);
thisMEGRDM = permute(groupAverage_MEG_RDM,[3 1 2]);
thisFMRIRDM = squeeze(FMRI_RDM(:,:,i_ROI));

%%

commonality2Models = zeros(numPicked2Models,lengthTrial);
sharedVariance1Model = zeros(numPicked3Models,lengthTrial);
RSquaredFmriMEG = zeros(1,lengthTrial);

permCommonality2Models = zeros([numPermutations,size(commonality2Models)]);
permSharedVariance1Model = zeros([numPermutations,size(sharedVariance1Model)]);
permRSquaredFmriMEG = zeros([numPermutations,size(RSquaredFmriMEG)]);

lowerTriangleIndices = find(tril(ones(widthRDM),-1));

%% Shared Variance (Model-MEG-fMRI)
for i_model = 1:numPicked3Models
    targetModelIndex = picked3Models(i_model);
    targetModel = models(targetModelIndex,:,:);
    validIndices = lowerTriangleIndices(find(~isnan(targetModel(lowerTriangleIndices))));
    [sharedVariance1Model(i_model,:), RSquaredFmriMEG] = combinedNModels(thisMEGRDM,thisFMRIRDM,targetModel,validIndices);
    for i_perm = 1:numPermutations
        disp([num2str(i_model) ' ' num2str(i_perm)]);
        rp = megPermutations(i_perm,:);
        permMEGRDM = thisMEGRDM(:,rp,rp);
        [permSharedVariance1Model(i_perm,i_model,:), permRSquaredFmriMEG(i_perm,:)] = combinedNModels(permMEGRDM,thisFMRIRDM,targetModel,validIndices);
    end
end


%% 2-model commonality

Mod1 = squeeze(models(1,:,:));
Mod2 = squeeze(models(2,:,:));
Mod3 = squeeze(models(3,:,:));

validIndices = lowerTriangleIndices(find(~isnan(Mod1(lowerTriangleIndices)) ...
                                       & ~isnan(Mod2(lowerTriangleIndices)) ...
                                       & ~isnan(Mod3(lowerTriangleIndices))));

for i_model = 1:numPicked2Models
    targetModelIndex = picked2Models(i_model);
    targetModel = squeeze(models(targetModelIndex,:,:));
    nonTargetModelsIndex = setdiff(picked2Models,targetModelIndex);
    nonTargetModels = models(nonTargetModelsIndex,:,:);
    commonality2Models(i_model,:) = commonalityAnalysisNModels(thisMEGRDM,thisFMRIRDM,targetModel,nonTargetModels,validIndices);
    for i_perm = 1:numPermutations
        disp([num2str(i_model) ' ' num2str(i_perm)]);
        rp = megPermutations(i_perm,:);
        permMEGRDM = thisMEGRDM(:,rp,rp);
        permCommonality2Models(i_perm,i_model,:) = commonalityAnalysisNModels(permMEGRDM,thisFMRIRDM,targetModel,nonTargetModels,validIndices);
    end
end

save(['Data/Processed/Figure6/Commonality/FMRI_MEG_Model_Commonality_' num2str(i_ROI) '.mat'],'commonality*','sharedVariance1Model','RSquaredFmriMEG', ...
    'permCommonality*','permSharedVariance1Model','permRSquaredFmriMEG');

