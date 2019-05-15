function plotFigure5D
%% Description: generates Figure5D 

% data = .mat file containing plotting variables outlined below
%
% 1) loads saved data files 
% 2) plots figure using:
%       - inter_h_[1-3] = vector of significance; 1 x timepoints for 3
%       comparison types:
%       {(Stimulus,Recognition),(Stimulus,PreGrey),(Recognition,PreGrey)}
%       - trace_mean = mean separability, timepoints x image type
%       - trace_ste = standard error of separability, timepoints x image type
%       - time = vector of timepoints, 1 x timepoints
%
%       ex) h_[1-3], significance for 3 image types (Predisamb, Postdisamb and Gray)
%           trace_mean, average RSA trace of 550 timepoints and 3 image types
%           trace_ste, standard error of RSA trace of 550 timepoints and 3 image types 
%           time, timepoints in 10ms steps for trial epoch (1x550)
%
% Dependencies: ~/boundedline


load('Data/Processed/Figure5/statistics_Figure5D.mat');

close all
figure
cols = 'rgb';
inter_hs = nan(3,550);
inter_hs(1,postStimulusTimeIndices) = inter_h_1;
inter_hs(2,postStimulusTimeIndices) = inter_h_2;
inter_hs(3,postStimulusTimeIndices) = inter_h_3;

for i = 1:3
    boundedline(time,meanTrace(i,:),steTrace(i,:),'alpha',cols(i));
    hold on
end

cols2 = 'kmc';
for i = 1:3
    plot(time,inter_hs(i,:) * 0.08 + (i*0.002),'color',cols2(i),'linewidth',4);
end
lg = {'Stimulus','Stimulus', ...
      'Recognition','Recognition', ...
      '(PreGrey)','(PreGrey)', ...
      'p<0.05,(Stimulus,Recognition)','p<0.05, (Stimulus,PreGrey)','p<0.05, (Recognition,PreGrey)'};
  
legend(lg,'location','eastoutside');
axis([-0.3 3 -0.01 0.1]);
xlabel('Time (s)')
ylabel('Separability: OffDiag - Diag');

set(gcf,'position',[251   323   942   468])
