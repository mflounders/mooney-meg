%% MEG Separability, diagonal analysis, with permutation tests
% This script is used in the manuscript to generate Figure 5, panel C
%
% 1) presets
% 2) initialize data
% 3) calculate dissimilarity for each condition
% 4) calculate original statistics, run permutations tests, threshold clusters
% 5) plot using plotFigure5C.m, saves plotting data
% (in Data/Processed/Figure5/statistics_Figure5C.mat)

% Estimated Duration: 90 minutes

clc
clearvars
close all
%% Parameters


numSubjects = 18;

lengthTrial = 550; %100 Hz sampling rate
time = -2.6:0.01:2.9-0.01;
postStimulusTimeIndices = 261:550; % post stimulus indices

% Images
numImages = 33;
PreImageIndices = 1:33;
PostImageIndices = 34:66;
GreyImageIndices = 67:99;

% Intra RDM comparison types
num_intraRDM = 3;
IntraRDMTypes = {'StimulusBased','RecognitionBased','PreGreyControl'};
IntraRDMLabels = {'(Pre,Post)','(Post,Grey)','(Pre,Grey)'};
IntraRDM.StimulusBased.triangle1 = PreImageIndices;
IntraRDM.StimulusBased.triangle2 = PostImageIndices;
IntraRDM.RecognitionBased.triangle1 = PostImageIndices;
IntraRDM.RecognitionBased.triangle2 = GreyImageIndices;
IntraRDM.PreGreyControl.triangle1 = PreImageIndices;
IntraRDM.PreGreyControl.triangle2 = GreyImageIndices;

upperTriangleIndices = find(triu(ones(numImages),1));

% Cluster defining statistics
cluster_def_p_thresh = 0.1;
plotting_p_thresh = 0.05;
nreps = 5000;

load('Data/Raw/MEG_Separability.mat');

%% Identify RDM diagonal for analysis
% Initialize
trace_mean = zeros(lengthTrial,num_intraRDM);   %mean for each RDM diagonal comparison
trace_ste = zeros(lengthTrial,num_intraRDM);    %standard error of mean for each RDM diagonal comparison
significantClusters = struct;                   %Stores significant clusters for each RDM diagonal comparison
allDiags = struct;

for i_intraRDM = 1:num_intraRDM
    
    this_intraRDMType = IntraRDMTypes{i_intraRDM};
    triangle1 = IntraRDM.(this_intraRDMType).triangle1;
    triangle2 = IntraRDM.(this_intraRDMType).triangle2;
    
    allDiag = zeros(lengthTrial,numSubjects,numImages);
    
    for i_sub = 1:numSubjects
        diags = zeros(lengthTrial,numImages);
        
        for i_time = 1:lengthTrial
            data1 = squeeze(MEG_Separability(triangle1,triangle2,i_time,i_sub));
            data1diag = diag(data1);
            
            % subset diagonal of desired RDM
            diags(i_time,:) = data1diag;
            
        end
        % collect between-condition diagonals for all subs
        allDiag(:,i_sub,:) = diags;
        
    end
    % store between-condition diagonals for all types
    allDiags.(this_intraRDMType).diag = allDiag;
    clear data1 data1z
    
end

%% Calculate mean of diagonals and original statistics

% Intialize
lengthPostStimulusTimeIndices = length(postStimulusTimeIndices);
inter_p_timecourse_1_orig = zeros(lengthPostStimulusTimeIndices,1);
inter_p_timecourse_2_orig = zeros(lengthPostStimulusTimeIndices,1);
inter_p_timecourse_3_orig = zeros(lengthPostStimulusTimeIndices,1);

inter_stats_timecourse_1_orig = zeros(lengthPostStimulusTimeIndices,1);
inter_stats_timecourse_2_orig = zeros(lengthPostStimulusTimeIndices,1);
inter_stats_timecourse_3_orig = zeros(lengthPostStimulusTimeIndices,1);

% Mean
x = mean(allDiags.StimulusBased.diag(postStimulusTimeIndices,:,:),3)';
y = mean(allDiags.RecognitionBased.diag(postStimulusTimeIndices,:,:),3)';
z = mean(allDiags.PreGreyControl.diag(postStimulusTimeIndices,:,:),3)';

% Statistics
for i_time = 1:lengthPostStimulusTimeIndices
    i_time
    % Stim, Recog
    [inter_p_timecourse_1_orig(i_time),~,stats] = signrank(x(:,i_time),y(:,i_time));
    inter_stats_timecourse_1_orig(i_time) = stats.signedrank - (sum(1:numSubjects)/2);
    
    % Stim,Pregrey
    [inter_p_timecourse_2_orig(i_time),~,stats2] = signrank(x(:,i_time),z(:,i_time));
    inter_stats_timecourse_2_orig(i_time) = stats2.signedrank  - (sum(1:numSubjects)/2);
    
    % Recog, PreGrey
    [inter_p_timecourse_3_orig(i_time),~,stats3] = signrank(y(:,i_time),z(:,i_time));
    inter_stats_timecourse_3_orig(i_time) = stats3.signedrank  - (sum(1:numSubjects)/2);
    
end

%% Statistics, pt.1 - inter-output, non-parametric comparisons 

% find significant clusters in original data time course
inter_clusters_1_orig = find_temporal_clusters(inter_stats_timecourse_1_orig, inter_p_timecourse_1_orig, cluster_def_p_thresh);
inter_clusters_2_orig = find_temporal_clusters(inter_stats_timecourse_2_orig, inter_p_timecourse_2_orig, cluster_def_p_thresh);
inter_clusters_3_orig = find_temporal_clusters(inter_stats_timecourse_3_orig, inter_p_timecourse_3_orig, cluster_def_p_thresh);
% permutation testing, finding distribution of maximum cluster size for shuffled data
[inter_clusters_shuffle_1, inter_shuffleMaxStat_1] = RSA_permutation_signrank(x,y,cluster_def_p_thresh,numSubjects,nreps);
[inter_clusters_shuffle_2, inter_shuffleMaxStat_2] = RSA_permutation_signrank(x,z,cluster_def_p_thresh,numSubjects,nreps);
[inter_clusters_shuffle_3, inter_shuffleMaxStat_3] = RSA_permutation_signrank(y,z,cluster_def_p_thresh,numSubjects,nreps);
% threshold significant clusters by comparing cluster statistics of original data with shuffled data
[inter_h_1] = temporal_cluster_thresh(inter_clusters_1_orig,inter_shuffleMaxStat_1,lengthPostStimulusTimeIndices,nreps);
[inter_h_2] = temporal_cluster_thresh(inter_clusters_2_orig,inter_shuffleMaxStat_2,lengthPostStimulusTimeIndices,nreps);
[inter_h_3] = temporal_cluster_thresh(inter_clusters_3_orig,inter_shuffleMaxStat_3,lengthPostStimulusTimeIndices,nreps);

%% Statistics, pt. 2 - output parametric testing against 0

% Initialize
p_timecourse_1_orig = zeros(lengthPostStimulusTimeIndices,1);
p_timecourse_2_orig = zeros(lengthPostStimulusTimeIndices,1);
p_timecourse_3_orig = zeros(lengthPostStimulusTimeIndices,1);

stats_timecourse_1_orig = zeros(lengthPostStimulusTimeIndices,1);
stats_timecourse_2_orig = zeros(lengthPostStimulusTimeIndices,1);
stats_timecourse_3_orig = zeros(lengthPostStimulusTimeIndices,1);

% find significant clusters in original data time course
for i_time = 1:lengthPostStimulusTimeIndices
    [~, p_timecourse_1_orig(i_time),~,statz] = ttest(x(:,i_time),0);
    stats_timecourse_1_orig(i_time) = statz.tstat;
    
    [~, p_timecourse_2_orig(i_time),~,statz] = ttest(y(:,i_time),0);
    stats_timecourse_2_orig(i_time) = statz.tstat;
    
    [~, p_timecourse_3_orig(i_time),~,statz] = ttest(z(:,i_time),0);
    stats_timecourse_3_orig(i_time) = statz.tstat;
end

clusters_orig_1 = find_temporal_clusters(stats_timecourse_1_orig, p_timecourse_1_orig, cluster_def_p_thresh);
clusters_orig_2 = find_temporal_clusters(stats_timecourse_2_orig, p_timecourse_2_orig, cluster_def_p_thresh);
clusters_orig_3 = find_temporal_clusters(stats_timecourse_3_orig, p_timecourse_3_orig, cluster_def_p_thresh);

% initialize for permutations
p_timecourse_1_shuf = zeros(nreps,lengthPostStimulusTimeIndices);
p_timecourse_2_shuf = zeros(nreps,lengthPostStimulusTimeIndices);
p_timecourse_3_shuf = zeros(nreps,lengthPostStimulusTimeIndices);
t_timecourse_1_shuf = zeros(nreps,lengthPostStimulusTimeIndices);
t_timecourse_2_shuf = zeros(nreps,lengthPostStimulusTimeIndices);
t_timecourse_3_shuf = zeros(nreps,lengthPostStimulusTimeIndices);
clusters_shuffle_1 = cell(nreps,1);
clusters_shuffle_2 = cell(nreps,1);
clusters_shuffle_3 = cell(nreps,1);
shuffleMaxStat_1 = zeros(nreps,1);
shuffleMaxStat_2 = zeros(nreps,1);
shuffleMaxStat_3 = zeros(nreps,1);
% permutation testing 
for i_rep = 1:nreps
    % compute a random vector for permuting subject labels
    random_vek = randi(2,1,numSubjects)';
    random_vek(random_vek==2) = -1;
    i_rep
    for i_time = 1:lengthPostStimulusTimeIndices
        % model data
        
        a = x(:,i_time);
        b = a.*random_vek;
        [~, p_timecourse_1_shuf(i_rep,i_time),~,statz] = ttest(b,0);
        t_timecourse_1_shuf(i_rep,i_time) = statz.tstat;  clear statz
        
        a = y(:,i_time);
        b = a.*random_vek;
        [~, p_timecourse_2_shuf(i_rep,i_time),~,statz] = ttest(b,0);
        t_timecourse_2_shuf(i_rep,i_time) = statz.tstat;  clear statz
        
        a = z(:,i_time);
        b = a.*random_vek;
        [~, p_timecourse_3_shuf(i_rep,i_time),~,statz] = ttest(b,0);
        t_timecourse_3_shuf(i_rep,i_time) = statz.tstat;  clear statz
        
    end
    
    % find distribution of maximum cluster size for shuffled data
    clusters_shuffle_1{i_rep} = find_temporal_clusters(t_timecourse_1_shuf(i_rep,:), p_timecourse_1_shuf(i_rep,:), cluster_def_p_thresh);
    shuffleMaxStat_1(i_rep) = clusters_shuffle_1{i_rep}.maxStatSumAbs; % two sided test
    
    clusters_shuffle_2{i_rep} = find_temporal_clusters(t_timecourse_2_shuf(i_rep,:), p_timecourse_2_shuf(i_rep,:), cluster_def_p_thresh);
    shuffleMaxStat_2(i_rep) = clusters_shuffle_2{i_rep}.maxStatSumAbs; % two sided test
    
    clusters_shuffle_3{i_rep} = find_temporal_clusters(t_timecourse_3_shuf(i_rep,:), p_timecourse_3_shuf(i_rep,:), cluster_def_p_thresh);
    shuffleMaxStat_3(i_rep) = clusters_shuffle_3{i_rep}.maxStatSumAbs; % two sided test
end
clear p_timecourse_shuf t_timecourse_shuf;

% Compare cluster statistics of original data with shuffled data
h_1 = nan(1,lengthPostStimulusTimeIndices);
for i_cluster = 1:clusters_orig_1.nClusters
    pval = sum(shuffleMaxStat_1 > clusters_orig_1.cluster_statSum(i_cluster)) / nreps;
    clusters_orig_1.cluster_pval(i_cluster) = pval;
    if clusters_orig_1.cluster_pval(1,i_cluster) < plotting_p_thresh
        timeCluster = clusters_orig_1.cluster_samples{i_cluster}(:,:);
        h_1(1,timeCluster) = 1;
    end
end

h_2 = nan(1,lengthPostStimulusTimeIndices);
for i_cluster = 1:clusters_orig_2.nClusters
    pval = sum(shuffleMaxStat_2 > clusters_orig_2.cluster_statSum(i_cluster)) / nreps;
    clusters_orig_2.cluster_pval(i_cluster) = pval;
    if clusters_orig_2.cluster_pval(1,i_cluster) < plotting_p_thresh
        timeCluster = clusters_orig_2.cluster_samples{i_cluster}(:,:);
        h_2(1,timeCluster) = 1;
    end
end

h_3 = nan(1,lengthPostStimulusTimeIndices);
for i_cluster = 1:clusters_orig_3.nClusters
    pval = sum(shuffleMaxStat_3 > clusters_orig_3.cluster_statSum(i_cluster)) / nreps;
    clusters_orig_3.cluster_pval(i_cluster) = pval;
    if clusters_orig_3.cluster_pval(1,i_cluster) < plotting_p_thresh
        timeCluster = clusters_orig_3.cluster_samples{i_cluster}(:,:);
        h_3(1,timeCluster) = 1;
    end
end

%% compute final diagonal means and ste for figure plot
meanTrace = zeros(num_intraRDM,lengthTrial);
steTrace = zeros(num_intraRDM,lengthTrial);
for i_intraRDM = 1:num_intraRDM
    this_intraRDMType = IntraRDMTypes{i_intraRDM};
    meanTrace(i_intraRDM,:) =  squeeze(mean(mean(allDiags.(this_intraRDMType).diag,3),2));
    x = mean((allDiags.(this_intraRDMType).diag),3);
    steTrace(i_intraRDM,:) =  (squeeze(std(x,[],2)) /  sqrt(numSubjects))';
end


save('Data/Processed/Figure5/statistics_Figure5C.mat','h_*','inter_h*','time','postStimulusTimeIndices','IntraRDMLabels','*Trace');

