function [clusters_shuffle, shuffleMaxStat] = RSA_permutation_signrank(data1,data2,p_thresh,nSubs,nReps)
% Description: cluster-based permutation testing of Wilcoxon signrank test
%
% input
% -----
% data1 = matrix of nSamples which holds dissimilarity values at each sample for a given image type, e.g. Pre-disambiguation images
% data2 = matrix of nSamples which holds dissimilarity values at each sample for a given image type, e.g. Post-disambiguation images
% p_thresh = p-value below which sensors will be considered for cluster analysis
% nsubs = number of subjects
% nreps = number of permutations for statistical testing

% output
% ------
% clusters_shuffle = struct included all statistical cluster information
% shuffleMaxStat = vector of nSamples which hold max cluster stat for thresholding and plotting significance at each sample

%% prep step
%remove singelton dimension
data1 = squeeze(data1);
data2 = squeeze(data2);
% pre allocate for speed
p_timecourse(nReps,size(data1,3)) = NaN;
stats_timecourse(nReps,size(data1,3)) = NaN;
clusters_shuffle{nReps} = NaN;
shuffleMaxStat(nReps) = NaN;

%% code

for i_rep = 1:nReps
    i_rep
shuffled_data_matrix(18,2) = NaN;
% compute a random vector for permuting subject labels
random_vek = randi(2,1,nSubs)';    
    
    for i_time = 1:size(data1,2) % should be same length as data2*
        clear h p ci stats a b
            data_matrix(:,1) = data1(:,i_time); 
            data_matrix(:,2) = data2(:,i_time);
            for i_sub = 1:length(data_matrix)
                if random_vek(i_sub) == 2
                    shuffled_data_matrix(i_sub,:) = fliplr(data_matrix(i_sub,:));
                else
                    shuffled_data_matrix(i_sub,:) = data_matrix(i_sub,:);
                end
            end 
        % paired signrank across subjects
        [p,h,stats] = signrank(shuffled_data_matrix(:,2),shuffled_data_matrix(:,1));
        p_timecourse(i_rep,i_time) = p;
        stats_timecourse(i_rep,i_time) = stats.signedrank - (sum(1:18)/2);
    end
    
    
    clusters_shuffle{i_rep} = find_temporal_clusters(stats_timecourse(i_rep,:), p_timecourse(i_rep,:), p_thresh);
    shuffleMaxStat(i_rep) = clusters_shuffle{i_rep}.maxStatSumAbs; % two sided test
end

end
