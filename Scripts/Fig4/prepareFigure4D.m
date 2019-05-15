%% MEG Representational Similarity Analysis (RSA) group plotting, with permutation tests
% This script is used in the manuscript to generate Figure 4, panel D 
%
% 1) Parameters
% 2) initialize data
% 3) Calculate mean dissimilarity for each subject and presentation type
% 4) calculate original statistics, run permutations tests, threshold clusters
% 5) saves data necessary to plot

% Estimated Duration: 180 minutes

clearvars
clc
%% Parameters
numSubjects = 18;
lengthTrial = 550;
time = -2.6:0.01:2.9-0.01;

%images
numImages = 33;
PreImageIndices = 1:33;
PostImageIndices = 34:66;
GreyImageIndices = 67:99;
num_presentationTypes = 3;
upperTriangleIndices = find(triu(ones(numImages),1));

%stats
nreps = 5000;
cluster_def_p_thresh = 0.1;

%plotting
save_path = 'Data/Processed/Figure4/';
save_header = 'statistics_Figure4D'; % saved as .mat file of variables

%% Initialize data
load('Data/Raw/MEG_RDM.mat');

%% Calculate mean dissimilarity for each subject and presentation type

dissim_trace = zeros(numSubjects,num_presentationTypes,lengthTrial);
for i_time = 1:lengthTrial
    for i_sub = 1:numSubjects
        Pre_RSA_matrix = squeeze(MEG_RDM(PreImageIndices,PreImageIndices,i_time,i_sub));
        Pre_mean = mean(Pre_RSA_matrix(upperTriangleIndices));
        
        Post_RSA_matrix = squeeze(MEG_RDM(PostImageIndices,PostImageIndices,i_time,i_sub));
        Post_mean = mean(Post_RSA_matrix(upperTriangleIndices));
        
        Grey_RSA_matrix = squeeze(MEG_RDM(GreyImageIndices,GreyImageIndices,i_time,i_sub));
        Grey_mean = mean(Grey_RSA_matrix(upperTriangleIndices));
        
        dissim_trace(i_sub,1,i_time) = Pre_mean;
        dissim_trace(i_sub,2,i_time) = Post_mean;
        dissim_trace(i_sub,3,i_time) = Grey_mean;
    end
end

trace_ste = zeros(lengthTrial,num_presentationTypes);
trace_mean = zeros(lengthTrial,num_presentationTypes);
for i_presentationType = 1:num_presentationTypes
    trace_ste(:,i_presentationType) = (squeeze(std(dissim_trace(:,i_presentationType,:))) /  sqrt(numSubjects))';
    trace_mean(:,i_presentationType) = squeeze(mean(dissim_trace(:,i_presentationType,:)));
end

%% calculate original statistics of dissimilarity across subjects

p_timecourse_1_orig = zeros(lengthTrial,1);
p_timecourse_2_orig = zeros(lengthTrial,1);
p_timecourse_3_orig = zeros(lengthTrial,1);

stats_timecourse_1_orig = zeros(lengthTrial,1);
stats_timecourse_2_orig = zeros(lengthTrial,1);
stats_timecourse_3_orig = zeros(lengthTrial,1);

for i_time = 1:lengthTrial
    
    % Post>Pre wilcoxon signed rank
    [p_timecourse_1_orig(i_time),~,stats] = signrank(dissim_trace(:,2,i_time),dissim_trace(:,1,i_time));
    stats_timecourse_1_orig(i_time) = stats.signedrank - (sum(1:numSubjects)/2);
    
    % Grey>Pre wilcoxon signed rank
    [p_timecourse_2_orig(i_time),~,stats2] = signrank(dissim_trace(:,3,i_time),dissim_trace(:,1,i_time));
    stats_timecourse_2_orig(i_time) = stats2.signedrank - (sum(1:numSubjects)/2);
    
    % Grey>Post wilcoxon signed rank
    [p_timecourse_3_orig(i_time),~,stats3] = signrank(dissim_trace(:,3,i_time),dissim_trace(:,2,i_time));
    stats_timecourse_3_orig(i_time) = stats3.signedrank - (sum(1:numSubjects)/2);
end

% find significant clusters in original data time course
clusters_1_orig = find_temporal_clusters(stats_timecourse_1_orig, p_timecourse_1_orig, cluster_def_p_thresh);
clusters_2_orig = find_temporal_clusters(stats_timecourse_2_orig, p_timecourse_2_orig, cluster_def_p_thresh);
clusters_3_orig = find_temporal_clusters(stats_timecourse_3_orig, p_timecourse_3_orig, cluster_def_p_thresh);


%% permutation testing

[clusters_shuffle_1, shuffleMaxStat_1] = RSA_permutation_signrank(dissim_trace(:,1,:),dissim_trace(:,2,:),cluster_def_p_thresh,numSubjects,nreps);
[clusters_shuffle_2, shuffleMaxStat_2] = RSA_permutation_signrank(dissim_trace(:,1,:),dissim_trace(:,3,:),cluster_def_p_thresh,numSubjects,nreps);
[clusters_shuffle_3, shuffleMaxStat_3] = RSA_permutation_signrank(dissim_trace(:,2,:),dissim_trace(:,3,:),cluster_def_p_thresh,numSubjects,nreps);


% test and calculate p-values against null distribution, plotting pval threshold =
% 0.05, which is set within the "temporal_cluster_thresh" function

[h_1] = temporal_cluster_thresh(clusters_1_orig,shuffleMaxStat_1,lengthTrial,nreps);
[h_2] = temporal_cluster_thresh(clusters_2_orig,shuffleMaxStat_2,lengthTrial,nreps);
[h_3] = temporal_cluster_thresh(clusters_3_orig,shuffleMaxStat_3,lengthTrial,nreps);




%% plotting

save(fullfile(save_path,save_header),'h_1', 'h_2', 'h_3', 'trace_mean', 'trace_ste', 'time');