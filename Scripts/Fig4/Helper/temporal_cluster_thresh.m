function [sig_vector] = temporal_cluster_thresh(clusters_orig,shuffleMaxStat,timewin,nReps)
% Description: threshold original stat clusters using maximum statistic from permuted null distribution 
%
% input
% -----
% clusters_orig = struct of original cluster statistics, nSamples in length
% shuffleMaxStat = vector of nSamples which hold max cluster stat for thresholding and plotting significance at each sample
% timewin = value of nSamples
% nReps = value of permutations 
%
% output
% ------
% sig_vector = vector of nSamples which hold a value of NaN or 1 if cluster is n.s. or significant, respectively 
%

%% code

sig_vector(1,timewin) = zeros;
for i_cluster = 1:clusters_orig.nClusters
    pval = sum(shuffleMaxStat > abs(clusters_orig.cluster_statSum(i_cluster)) ) / nReps;
    clusters_orig.cluster_pval(i_cluster) = pval;
    if clusters_orig.cluster_pval(1,i_cluster) < 0.05
        highlite_sens = clusters_orig.cluster_samples{i_cluster}(:,:);
        sig_vector(1,highlite_sens) = 1;
    end
    sig_vector(sig_vector==0)=nan;
end

end
