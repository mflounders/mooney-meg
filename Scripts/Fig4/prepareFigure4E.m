%% MEG Representational Similarity Analysis (RSA) intra-RDM analysis, with permutation tests
% This script is used in the manuscript to generate Figure 4, panel E
%
% 1) presets
% 2) initialize data
% 3) calculate dissimilarity for each condition
% 4) calculate original statistics, run permutations tests, threshold clusters
% 5) saves data necessary to plot
% (in Data/Processed/Figure4/statistics_Figure4E.mat)

% Estimated Duration: 40-60 minutes

clc
clearvars

%% parameters

numSubjects = 18;
lengthTrial = 550; %100 Hz sampling rate 
time = -2.6:0.01:2.9-0.01;

% Images
num_images = 33;
PreImageIndices = 1:33;
PostImageIndices = 34:66;
GreyImageIndices = 67:99;

% Intra RDM comparison types
numIntraRDM = 2;
IntraRDMTypes = {'StimulusBased','RecognitionBased'};
IntraRDM.StimulusBased.triangle1 = PreImageIndices;
IntraRDM.StimulusBased.triangle2 = PostImageIndices;
IntraRDM.RecognitionBased.triangle1 = PostImageIndices;
IntraRDM.RecognitionBased.triangle2 = GreyImageIndices;

upperTriangleIndices = find(triu(ones(num_images),1));

% Cluster defining statistics
cluster_def_p_thresh = 0.1;
plotting_p_thresh = 0.05;
nreps = 5000;

%plotting
save_path = 'Data/Processed/Figure4/';
save_header = 'statistics_Figure4E'; % saved as .mat file of variables

load('Data/Raw/MEG_RDM.mat');

%% Initialize
trace_mean = zeros(lengthTrial,numIntraRDM);   %mean for each intra_RDM comparison
trace_ste = zeros(lengthTrial,numIntraRDM);    %standard error of mean for each intra_RDM comparison
significantClusters = struct;                   %Stores significant clusters for each intra_RDM comparison

for i_intraRDM = 1:numIntraRDM
    
    this_intraRDMType = IntraRDMTypes{i_intraRDM};
    triangle1 = IntraRDM.(this_intraRDMType).triangle1;
    triangle2 = IntraRDM.(this_intraRDMType).triangle2;
    
    %% Calculate mean and sem for this intra-rdm analysis 
    pearz = zeros(lengthTrial,numSubjects);
    for i_sub = 1:numSubjects
        rho = zeros(lengthTrial,1);
        for i_time = 1:lengthTrial
            data1 = 1-squeeze(MEG_RDM(triangle1,triangle1,i_time,i_sub));   %calculate Similarity from dissimilarity
            data1z = pear_fisherz(data1(upperTriangleIndices));
            
            data2 = 1-squeeze(MEG_RDM(triangle2,triangle2,i_time,i_sub));
            data2z = pear_fisherz(data2(upperTriangleIndices));
            
            % Pearson correlation between model sections
            rho(i_time) = corr(data1z,data2z,'Type','Pearson');
        end
        % transform Pearson's R to Z
        pearz(:,i_sub) = pear_fisherz(rho');
    end
    trace_ste(:,i_intraRDM)= (squeeze(std(pearz,[],2)) /  sqrt(numSubjects))';
    trace_mean(:,i_intraRDM)= mean(pearz,2);
    
    clear data1 data1z data2 data2z rho
       
    %% Find significant clusters for unshuffled data
    p_timecourse_orig = zeros(lengthTrial,1);
    stats_timecourse_orig = zeros(lengthTrial,1);
    for i_time = 1:lengthTrial
        % one sample ttest against 0
        [~, p_timecourse_orig(i_time),~,statz] = ttest(pearz(i_time,:)',0);
        stats_timecourse_orig(i_time) = statz.tstat;
    end
    
    clusters_orig = find_temporal_clusters(stats_timecourse_orig, p_timecourse_orig, cluster_def_p_thresh);
    
    clear statz p_timecourse_orig stats_timecourse_orig;
    
    %% Find distribution of maximum cluster size for shuffled data
    
    p_timecourse_shuf = zeros(nreps,lengthTrial);
    t_timecourse_shuf = zeros(nreps,lengthTrial);
    clusters_shuffle = cell(nreps,1);
    shuffleMaxStat = zeros(nreps,1);
    for i_rep = 1:nreps
        % compute a random vector for permuting subject labels
        random_vek = randi(2,1,numSubjects)';
        random_vek(random_vek==2) = -1;
        
        for i_time = 1:lengthTrial
            % model data
            a = pearz(i_time,:)';
            % data with random vector of 1s and -1s applied to it
            b = a.*random_vek;
            
            % stat timecourse, one sample ttest against 0
            [~, p_timecourse_shuf(i_rep,i_time),~,statz] = ttest(b,0);
            t_timecourse_shuf(i_rep,i_time) = statz.tstat;  clear statz
        end
        
        clusters_shuffle{i_rep} = find_temporal_clusters(t_timecourse_shuf(i_rep,:), p_timecourse_shuf(i_rep,:), cluster_def_p_thresh);
        shuffleMaxStat(i_rep) = clusters_shuffle{i_rep}.maxStatSumAbs; % two sided test
    end
    clear p_timecourse_shuf t_timecourse_shuf;
    
    %% Compare cluster statistics of original data, to shuffled data
    significantClusters.(this_intraRDMType).timeCourse = nan(1,lengthTrial);
    for i_cluster = 1:clusters_orig.nClusters
        pval = sum(shuffleMaxStat > clusters_orig.cluster_statSum(i_cluster)) / nreps;
        clusters_orig.cluster_pval(i_cluster) = pval;
        if clusters_orig.cluster_pval(1,i_cluster) < plotting_p_thresh
            timeCluster = clusters_orig.cluster_samples{i_cluster}(:,:);
            significantClusters.(this_intraRDMType).timeCourse(1,timeCluster) = 1;
        end
    end
end
%% manuscript plot

save(fullfile(save_path,save_header),'significantClusters', 'trace_mean', 'trace_ste', 'time');



