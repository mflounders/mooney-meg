function MEG_Model_Fusion
%Correlates MEG RDMs with the three models.
%Also correlates row-shuffled MEG RDMs with the three models

%% PRESETS
lengthTrial = 550;
numTasks = 3;
numImages = 33;
numPermutations = 1000;
numModels = 3;
widthRDM = 99;
%%
load('Data/Raw/MEG_RDM.mat');
load('Data/Raw/MEG_Permutations.mat');   %Contains shuffled indices for each permutation
load(['Data/Raw/Models.mat']);   %Contains shuffled indices for each permutation

groupAverage_MEG_RDM = mean(MEG_RDM,4);

%%
thisMEGRDM = permute(groupAverage_MEG_RDM,[3 1 2]);
RModelMEG = zeros(numModels,lengthTrial); %Correlation of MEG RDM with each of the picked models
permRModelMEG = zeros([numPermutations,size(RModelMEG)]); %Correlation of shuffled MEG RDM with each of the picked models
lowerTriangleIndices = find(tril(ones(widthRDM),-1));


for i_model = 1:numModels
    targetModel = squeeze(models(i_model,:,:));
    notNans = find(~isnan(targetModel(lowerTriangleIndices)));
    for i_time = 1:lengthTrial
        timeMEGRDM = squeeze(thisMEGRDM(i_time,:,:));
        RModelMEG(i_model,i_time) = corr(timeMEGRDM(lowerTriangleIndices(notNans)),targetModel(lowerTriangleIndices(notNans)),'type','spearman');
    end
    
    for i_perm = 1:numPermutations
        disp([num2str(i_model) ' ' num2str(i_perm)]);
        rp = megPermutations(i_perm,:);
        permMEGRDM = thisMEGRDM(:,rp,rp);
        for i_time = 1:lengthTrial
            ptimeMEGRDM = squeeze(permMEGRDM(i_time,:,:));
            permRModelMEG(i_perm,i_model,i_time) = corr(ptimeMEGRDM(lowerTriangleIndices(notNans)),targetModel(lowerTriangleIndices(notNans)),'type','spearman');
        end
    end
end

save(['Data/Processed/Figure6/MEG_Model_Fusion.mat'],'RModelMEG','permRModelMEG');

