function plotFigure4D
%% Description: generates Figure4D

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


%% load up saved data and stats

load('Data/Processed/Figure4/statistics_Figure4D.mat');
disp('Data loaded, plotting');

%% plot Figure 4D

figure;
hold all
plt1 = boundedline(time,trace_mean(:,1),trace_ste(:,1),'alpha');
plt2 = boundedline(time,trace_mean(:,2),trace_ste(:,2),'alpha');
plt3 = boundedline(time,trace_mean(:,3),trace_ste(:,3),'alpha');

set(plt1,'Color',[0.815686285495758 0.670588254928589 0.325490206480026],'DisplayName','Pre','LineWidth',2.0);
set(plt2,'Color',[0 0.600000023841858 0.600000023841858],'DisplayName','Post','LineWidth',2.0);
set(plt3,'Color',[0 0 0],'DisplayName','Grey','LineWidth',1.5,'LineStyle',':');

plot(time,(h_1 + 0.01),'-*','DisplayName','Post trace, Pre trace,signrank');
plot(time,(h_2 + 0.03),'-*','DisplayName','Grey trace, Pre trace,signrank');
plot(time,h_3,'-*','DisplayName','Grey trace > Post trace,signrank');

ylabel('dissimilarity (1-r)','FontSize',14);
xlabel('time (s)','FontSize',14);
title('average RDM matrix traces','FontSize',14);
legend('show','location','northeastoutside');
xlim([-0.5 3]);
end