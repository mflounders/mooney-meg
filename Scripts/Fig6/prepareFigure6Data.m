function prepareFigure6Data
%Creates files for plotting Figures.  
%These files are already provided, so not necessary to run in order to plot
%figures
%
% Estimated duration: With 1000 permutations, this will take multiple days to complete

FMRI_Model_Fusion; %Correlates fMRI RDMs with the three models. Creates FMRI_Model_Fusion.mat
MEG_Model_Fusion;  %Correlates MEG RDMs and shuffled MEG RDMs with the three models. Creates MEG_Model_Fusion.mat
statistics_MEG_Model_Fusion;  %Calculates significant clusters for MEG-Model RDM. Creates statistics_Figure6D.mat

numROIs = 20;
for i = 1:numROIs
    FMRI_MEG_Model_Commonality(i_ROI);    %Calculates commonality analysis between MEG, fMRI, and model RDMs (+ permutations). Creates 'FMRI_MEG_Model_Commonality_[i_ROI].mat'
    statistics_FMRI_MEG_Model_Commonality(i_ROI);      %Calculates significant clusters for commonality analysis. Creates 'statistics_FMRI_MEG_Model_Commonality_[i_ROI].mat
end