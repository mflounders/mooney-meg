function statistics_FMRI_MEG_Model_Commonality(i_ROI)
%Calculate Significant clusters of commonality

%% Parameters
numModels = 3;
percentile95 = 0.05;
postStimulusTimes = 261:550;

load(['Data/Processed/Figure6/Commonality/FMRI_MEG_Model_Commonality_' num2str(i_ROI) '.mat'])
load('Data/Raw/FMRI_RDM.mat','ROI_labels');
numPermutations = size(permCommonality2Models,1);
numTimeSteps = size(commonality2Models,2);
cutoffCommonality2Models = zeros(numTimeSteps,2);
meanCommonality2Models = zeros(numTimeSteps,2);
cutoffSharedVariance1Model = zeros(numTimeSteps,3);
meanSharedVariance1Model = zeros(numTimeSteps,3);

%% Get cutoff and variance from permutations
tmp = squeeze(permRSquaredFmriMEG);
tmp_sorted = sort(tmp,'descend')';
cutoffRSquaredFmriMEG = tmp_sorted(:,floor(percentile95*size(tmp,1)));
meanRSquaredFmriMEG = tmp_sorted(:,floor(0.5*size(tmp,1)));
for i_model = 1:3    
    tmp = squeeze(permSharedVariance1Model(:,i_model,:));
    tmp_sorted = sort(tmp,'descend')';
    cutoffSharedVariance1Model(:,i_model) = tmp_sorted(:,floor(percentile95*size(tmp,1)));
    meanSharedVariance1Model(:,i_model) = tmp_sorted(:,floor(0.5*size(tmp,1)));
end
for i_model = 1:2
    tmp = squeeze(permCommonality2Models(:,i_model,:));
    tmp_sorted = sort(tmp,'descend')'; %#ok<*UDIM>
    cutoffCommonality2Models(:,i_model) = tmp_sorted(:,floor(percentile95*size(tmp,1)));
    meanCommonality2Models(:,i_model) = tmp_sorted(:,floor(0.5*size(tmp,1)));
    
end

clustersPermCommonality2Models = zeros(numPermutations,2);
clustersPermSharedVariance1Model = zeros(numPermutations,numModels);
clustersRSquaredFmriMEG = zeros(numPermutations,1);

for i_perm = 1:numPermutations
    % find clusters gives us the cluster sizes of the current permutation
    % in the current ROI that are larger than the cutoff value
    for i_model = 1:numModels
        clustersPermSharedVariance1Model(i_perm,i_model) = max(find_clusters(squeeze(permSharedVariance1Model(i_perm,i_model,postStimulusTimes))>cutoffSharedVariance1Model(postStimulusTimes,i_model)));
    end
    for i_model = 1:2
        clustersPermCommonality2Models(i_perm,i_model) = max(find_clusters(squeeze(permCommonality2Models(i_perm,i_model,postStimulusTimes))>cutoffCommonality2Models(postStimulusTimes,i_model)));
    end
    clustersRSquaredFmriMEG(i_perm) = max(find_clusters(squeeze(permRSquaredFmriMEG(i_perm,:,postStimulusTimes)) >cutoffRSquaredFmriMEG(postStimulusTimes)));
end

c_sorted = sort(clustersRSquaredFmriMEG(:),'descend');
clusterCutoffRSquaredFmriMEG = c_sorted(floor(percentile95*numPermutations));

for i_model = 1:numModels
    c_sorted = sort(clustersPermSharedVariance1Model(:,i_model),'descend');
    clusterCutoffSharedVariance1Model(i_model) = c_sorted(floor(percentile95*numPermutations));
end

for i_model = 1:2
    c_sorted = sort(clustersPermCommonality2Models(:,i_model),'descend');
    clusterCutoffCommonality2Models(i_model) = c_sorted(floor(percentile95*numPermutations));
end

clusterTimeSeriesCommonality2Models = zeros(numTimeSteps,2);
clusterTimeSeriesSharedVariance1Model = zeros(numTimeSteps,numModels);
clusterTimeSeriesRSquaredFmriMEG = zeros(numTimeSteps,1);

[c,~,~,clustind] = find_clusters(RSquaredFmriMEG(postStimulusTimes)>cutoffRSquaredFmriMEG(postStimulusTimes));
for i_c = 1:length(c) % loop over cluster sizes
    if c(i_c)<=clusterCutoffRSquaredFmriMEG% was <clust_cutoff_task(j_roi)
        clustind(clustind==i_c) = 0;
    end
end

clusterTimeSeriesRSquaredFmriMEG(postStimulusTimes) = clustind;
for i_model = 1:3    
    [c,~,~,clustind] = find_clusters(sharedVariance1Model(i_model,postStimulusTimes)'>cutoffSharedVariance1Model(postStimulusTimes,i_model));
    for i_c = 1:length(c) % loop over cluster sizes
        if c(i_c)<=clusterCutoffSharedVariance1Model(i_model)% was <clust_cutoff_task(j_roi)
            clustind(clustind==i_c) = 0;
        end
    end
    clusterTimeSeriesSharedVariance1Model(postStimulusTimes,i_model) = clustind;
end


for i_model = 1:2
    [c,~,~,clustind] = find_clusters(commonality2Models(i_model,postStimulusTimes)'>cutoffCommonality2Models(postStimulusTimes,i_model));
    for i_c = 1:length(c) % loop over cluster sizes
        if c(i_c)<=clusterCutoffCommonality2Models(i_model)% was <clust_cutoff_task(j_roi)
            clustind(clustind==i_c) = 0;
        end
    end
    clusterTimeSeriesCommonality2Models(postStimulusTimes,i_model) = clustind;
    
end
save(['Data/Processed/Figure6/Commonality/statistics_FMRI_MEG_Model_Commonality_' num2str(i_ROI) '.mat'],'clusterTimeSeries*','commonality*','sharedVariance1Model','RSquaredFmriMEG','ROI_labels');
