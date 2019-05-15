function plotFigure4E
%% Description: generates Figure4E 

% data = .mat file containing 4 plotting variables outlined below
%
% 1) loads saved data files 
% 2) plots figure using:
%       - h_[imagetype] = vector of significance; 1 x timepoints 
%       - trace_mean = mean representation dissimilarity, timepoints x image type
%       - trace_ste = standard error of dissimilarity, timepoints x image type
%       - time = vector of timepoints, 1 x timepoints
%
%       ex) h_[1-3], significance for 3 image types (Predisamb, Postdisamb and Gray)
%           trace_mean, average RSA trace of 550 timepoints and 3 image types
%           trace_ste, standard error of RSA trace of 550 timepoints and 3 image types 
%           time, timepoints in 10ms steps for trial epoch (1x550)
%
% Dependencies: ~/boundedline

%% load data

    load('Data/Processed/Figure4/statistics_Figure4E.mat');
    disp('Data loaded, plotting saved files');


%% plot

color1 = 'g';
color2 = 'm';
figure;
plt1 = boundedline(time,trace_mean(:,1),trace_ste(:,1),'alpha');
hold all
plt2 = boundedline(time,trace_mean(:,2),trace_ste(:,2),'alpha');
set(plt1,'Color',color1,'DisplayName','Stimulus-based','LineWidth',2.0);
set(plt2,'Color',color2,'DisplayName','Recognition-based','LineWidth',2.0);
legend('show','location','northeastoutside');
plot(time,(significantClusters.StimulusBased.timeCourse * 0.19),'-*','DisplayName','Stim-Based, ttest','color',color1);
plot(time,(significantClusters.RecognitionBased.timeCourse * 0.18),'-*','DisplayName','Recog-Based, ttest','color',color2);

xlim([-0.5 3]);
ylabel('Fisher Z','FontSize',14);
xlabel('time (s)','FontSize',14);
title('intra-RDM analysis','FontSize',14);