function plotFigure5C
%% Description: generates Figure5C

% data = .mat file containing plotting variables outlined below
%
% 1) loads saved data files 
% 2) plots figure using:
%       - h_[1-3], significance for (Stimulus, Post and PreGrey): 1 x timepoints  
%       - inter_h_[1-3] = vector of significance; 1 x timepoints for 3
%       comparison types:
%       {(Stimulus,Recognition),(Stimulus,PreGrey),(Recognition,PreGrey)}
%       - trace_mean = mean separability, timepoints x image type
%       - trace_ste = standard error of separability, timepoints x image type
%       - time = vector of timepoints, 1 x timepoints
%
% Dependencies: ~/boundedline

load('Data/Processed/Figure5/statistics_Figure5C.mat');

%% Figure Parameters
cols = 'RGB';
cols2 = 'kmc';
st_h = 0.08;
st_h2 = 0.09;
d_h = 0.003;
ylim = [-0.02 0.11];
lg = {'Stimulus','Stimulus','p<0.05, Stimulus','Recognition','Recognition','p<0.05, Recognition','PreGrey','PreGrey','p<0.05, PreGrey'};
lg{10} = 'p<0.05, (Stimulus,Recognition)';
lg{11} = 'p<0.05, (Stimulus,PreGrey)';
lg{12} = 'p<0.05, (Recognition,PreGrey)';

%%

hs = nan(3,550);
inter_hs = nan(3,550);

hs(1,postStimulusTimeIndices) = h_1;
hs(2,postStimulusTimeIndices) = h_2;
hs(3,postStimulusTimeIndices) = h_3;
inter_hs(1,postStimulusTimeIndices) = inter_h_1;
inter_hs(2,postStimulusTimeIndices) = inter_h_2;
inter_hs(3,postStimulusTimeIndices) = inter_h_3;

time = -2.6:0.01:2.9-0.01;
%%
figure

for i_intraRDM = 1:3
    boundedline(time,meanTrace(i_intraRDM,:),steTrace(i_intraRDM,:),'alpha',cols(i_intraRDM));
    hold on
    plot(time,hs(i_intraRDM,:) * st_h2 + (i_intraRDM*d_h),'color',cols(i_intraRDM),'linewidth',4);    
end

for i = 1:3
    plot(time,inter_hs(i,:) * st_h + (i*d_h),'color',cols2(i),'linewidth',4);
end

legend(lg,'location','EastOutside')
xlabel('Time (seconds)')
ylabel('rz(within)-rz(between)');
title('rz(within)-rz(between)')
axis([-0.5 3 ylim(1) ylim(2)])
set(gcf,'position',[251   323   942   468])
