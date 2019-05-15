%% commonality analysis plotting - Figure 6E 
% significance plots for all ROIs

numROIs = 20;
lengthTrial = 550;
numModels = 3;
plotTimePostStim = 3;
tms = -2.6:0.01:2.9-0.01;   %time (seconds) since stimulus 
c2Models = [1 3];
comms2Model  = zeros(numModels,lengthTrial,2);
sharedVariance = zeros(20,lengthTrial,numModels);
load(['Data/Processed/Figure6/statistics_Figure6D.mat']);
load(['Data/Processed/Figure6/FRMI_Model_Fusion.mat']);
%% Load Data
for i_ROI = 1:numROIs
    load(['Data/Processed/Figure6/Commonality/statistics_FMRI_MEG_Model_Commonality_' num2str(i_ROI) '.mat']);
    
    for i_model = 1:2
        comms2Model(i_ROI,:,i_model) = clusterTimeSeriesCommonality2Models(:,i_model) > 0 & RModelMEG(c2Models(i_model),:)'>0 & rFMRIModel(c2Models(i_model),i_ROI)'>0;
    end
    
    for i_model = 1:3
        sharedVariance(i_ROI,:,i_model) = clusterTimeSeriesSharedVariance1Model(:,i_model) > 0 & RModelMEG(i_model,:)'>0 & rFMRIModel(i_model,i_ROI)'>0;
    end
end

legendText = ROI_labels;


%% Plot Figures
figure
    modelNames = {'Stimulus','Attention','Recognition'};
    imagesc(tms,1:numROIs,squeeze(comms2Model(:,:,1))>0,[0 1]);
    axis([0 plotTimePostStim 0 numROIs]);
    set(gca,'ytick',[1:numROIs])
    set(gca,'yticklabel',legendText)
    title(['Commonality 2Model ' modelNames{1}]);    
    
figure;
    imagesc(tms,1:numROIs,squeeze(comms2Model(:,:,2))>0,[0 1]);
    axis([0 plotTimePostStim 0 numROIs]);
    set(gca,'ytick',[1:numROIs])
    set(gca,'yticklabel',legendText)
    title(['Commonality 2Model ' modelNames{3}]);
    
figure;
    imagesc(tms,1:numROIs,squeeze(sharedVariance(:,:,2))>0,[0 1]);
    axis([0 plotTimePostStim 0 numROIs]);
    set(gca,'ytick',[1:numROIs])
    set(gca,'yticklabel',legendText)
    xlabel('Time (seconds))');
    title(['SharedVariance ' modelNames{2}]);
