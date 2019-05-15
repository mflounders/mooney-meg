function statistics_MEG_Model_Fusion
%Calculates significant clusters for MEG_Model_Fusion

%% Parameters
lengthTrial = 550;
postStimulusTimes = 261:550;
numModels = 3;
percentile95 = 0.05;

load(['Data/Processed/Figure6/MEG_Model_Fusion.mat']);
numPermutations = size(permRModelMEG,1);

%% Find 95 percentile of cluster size for each model
cutoffRModelMEG = zeros(lengthTrial,numModels);
meanRModelMEG = zeros(lengthTrial,numModels);
for i_model = 1:numModels
    tmp = squeeze(permRModelMEG(:,i_model,:));
    tmp_sorted = sort(tmp,'descend')'; 
    cutoffRModelMEG(:,i_model) = tmp_sorted(:,floor(percentile95*size(tmp,1)));
    meanRModelMEG(:,i_model) = tmp_sorted(:,floor(0.5*size(tmp,1)));
end

clustersPermRModelMEG = zeros(numPermutations,numModels);
for i_perm = 1:numPermutations
    % find clusters gives us the cluster sizes of the current permutation
    % in the current ROI that are larger than the cutoff value
    for i_model = 1:numModels
        clustersPermRModelMEG(i_perm,i_model) = max(find_clusters(squeeze(permRModelMEG(i_perm,i_model,postStimulusTimes))>cutoffRModelMEG(postStimulusTimes,i_model)));    
    end
end

clusterCutoffRModelMEG = zeros(1,numModels);
for i_model = 1:numModels
    c_sorted = sort(clustersPermRModelMEG(:,i_model),'descend');
    clusterCutoffRModelMEG(i_model) = c_sorted(floor(percentile95*numPermutations));
end

%%
clusterTimeSeriesRModelMEG = zeros(lengthTrial,numModels);
for i_model = 1:numModels
    [c,~,~,clustind] = find_clusters(RModelMEG(i_model,postStimulusTimes)'>cutoffRModelMEG(postStimulusTimes,i_model));
    for i_c = 1:length(c) % loop over cluster sizes
        if c(i_c)<=clusterCutoffRModelMEG(i_model)% was <clust_cutoff_task(j_roi)
            clustind(clustind==i_c) = 0;
        end
    end
    clusterTimeSeriesRModelMEG(postStimulusTimes,i_model) = clustind;
end

save(['Data/Processed/statistics_Figure6D.mat'],'clusterTimeSeries*','RModelMEG');
