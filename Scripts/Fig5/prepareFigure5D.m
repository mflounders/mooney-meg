%% MEG Separability, off-diagonal vs. diagonal "strength" analysis, with permutation tests
% This script is used in the manuscript to generate Figure 5, panel D
%
% 1) presets
% 2) initialize data
% 3) calculate dissimilarity for each condition
% 4) calculate original statistics, run permutations tests, threshold clusters
% (in Data/Processed/Figure5/statistics_Figure5D.mat)

% Estimated Duration: 90 minutes

clc
clearvars
close all
%% Parameters

numSubjects = 18;
lengthTrial = 550; %100 Hz sampling rate
time = -2.6:0.01:2.9-0.01;
postStimulusTimeIndices = 261:550; %post stimulus indices

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
lowerTriangleIndices = find(tril(ones(numImages),-1));
offDiagIndices = [upperTriangleIndices; lowerTriangleIndices];
numOffDiagElements = length(offDiagIndices);
% Cluster defining statistics
cluster_def_p_thresh = 0.1;
plotting_p_thresh = 0.05;
nreps = 5000;

load('Data/Raw/MEG_Separability.mat');

%% Identify RDM diagonals and off-diagonals for analysis
% Initialize
trace_mean = zeros(lengthTrial,num_intraRDM);   %mean for each RDM comparison
trace_ste = zeros(lengthTrial,num_intraRDM);    %standard error of mean for each RDM comparison
significantClusters = struct;                   %Stores significant clusters for each RDM comparison
allDiags = struct;

for i_intraRDM = 1:num_intraRDM
    
    this_intraRDMType = IntraRDMTypes{i_intraRDM};
    triangle1 = IntraRDM.(this_intraRDMType).triangle1;
    triangle2 = IntraRDM.(this_intraRDMType).triangle2;
    
    allDiag = zeros(lengthTrial,numSubjects,numImages);
    allOffDiag = zeros(lengthTrial,numSubjects,numOffDiagElements);
    for i_sub = 1:numSubjects
        diags = zeros(lengthTrial,numImages);
        offDiags = zeros(lengthTrial,numOffDiagElements);
        for i_time = 1:lengthTrial
            data1 = squeeze(MEG_Separability(triangle1,triangle2,i_time,i_sub));   %calculate Similarity from dissimilarity
            data1diag = diag(data1);
            data1Offdiag = data1(offDiagIndices);
            % subset between-condition diagonal and off diagonal of desired RDM
            diags(i_time,:) = data1diag;
            offDiags(i_time,:) = data1Offdiag;
        end
        % collect between-condition diagonals and off diagonals for all subs
        allDiag(:,i_sub,:) = diags;
        allOffDiag(:,i_sub,:) = offDiags;
    end
    % store between-condition diagonals and off diagonals for all types
    allDiags.(this_intraRDMType).diag = allDiag;
    allDiags.(this_intraRDMType).offDiag = allOffDiag;
    
    clear data1 data1z
    
end

%% Calculate difference of means of off-diagonals and diagonals and original statistics
% Mean
IntraRDMTypes = {'StimulusBased','RecognitionBased','PreGreyControl'};
x = mean(allDiags.StimulusBased.offDiag,3) - mean(allDiags.StimulusBased.diag,3);
y = mean(allDiags.RecognitionBased.offDiag,3) - mean(allDiags.RecognitionBased.diag,3);
z = mean(allDiags.PreGreyControl.offDiag,3) - mean(allDiags.PreGreyControl.diag,3);

x = x(postStimulusTimeIndices,:)';
y = y(postStimulusTimeIndices,:)';
z = z(postStimulusTimeIndices,:)';

% Initialize
lengthPostStimulusTimeIndices = length(postStimulusTimeIndices);
inter_p_timecourse_1_orig = zeros(lengthPostStimulusTimeIndices,1);
inter_p_timecourse_2_orig = zeros(lengthPostStimulusTimeIndices,1);
inter_p_timecourse_3_orig = zeros(lengthPostStimulusTimeIndices,1);

inter_stats_timecourse_1_orig = zeros(lengthPostStimulusTimeIndices,1);
inter_stats_timecourse_2_orig = zeros(lengthPostStimulusTimeIndices,1);
inter_stats_timecourse_3_orig = zeros(lengthPostStimulusTimeIndices,1);

% Statistics
for i_time = 1:lengthPostStimulusTimeIndices
    i_time
    % Stim, Recog
    [inter_p_timecourse_1_orig(i_time),~,stats] = signrank(x(:,i_time),y(:,i_time));
    inter_stats_timecourse_1_orig(i_time) = stats.signedrank - (sum(1:numSubjects)/2);
    
    % Stim,Pregrey
    [inter_p_timecourse_2_orig(i_time),~,stats2] = signrank(x(:,i_time),z(:,i_time));
    inter_stats_timecourse_2_orig(i_time) = stats2.signedrank - (sum(1:numSubjects)/2);
    
    % Recog, PreGrey
    [inter_p_timecourse_3_orig(i_time),~,stats3] = signrank(y(:,i_time),z(:,i_time));
    inter_stats_timecourse_3_orig(i_time) = stats3.signedrank - (sum(1:numSubjects)/2);
    
end

%% Statistics - inter-output, non-parametric comparisons

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

%% compute final means and ste for figure plot
meanTrace = zeros(num_intraRDM,lengthTrial);
steTrace = zeros(num_intraRDM,lengthTrial);
for i_intraRDM = 1:num_intraRDM
    x = mean(allDiags.(IntraRDMTypes{i_intraRDM}).offDiag,3) - mean(allDiags.(IntraRDMTypes{i_intraRDM}).diag,3);
    meanTrace(i_intraRDM,:) =  mean(x,2);
    steTrace(i_intraRDM,:) =  (squeeze(std(x,[],2)) /  sqrt(numSubjects))';
    [~,p(i_intraRDM,:)] = ttest(x');
end

save('Data/Processed/Figure5/statistics_Figure5D.mat','time','postStimulusTimeIndices','IntraRDMLabels','*Trace','inter*');
