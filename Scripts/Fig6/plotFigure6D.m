%% Model MEG Fusion analysis plotting - Figure 6D

load(['Data/Processed/Figure6/statistics_Figure6D.mat']);

%close all
colors = 'rbyk';
tms = -2.6:0.01:2.9-0.01;
figure
set(gcf,'position',[37         127        1416/2         718])
for i_model = 1:3
    plot(tms,RModelMEG(i_model,:),'color',colors(i_model));
    hold on
    pcks = find(clusterTimeSeriesRModelMEG(:,i_model));
    if ~isempty(pcks)
        plot(tms(pcks), (0.4-(i_model*0.01)) * ones(length(pcks),1),'.','color',colors(i_model));
    else
        plot(0, 0,'.');
    end
end
legend({'StimulusModel-MEG','StimulusModel-MEG',...
    'AttentionModel-MEG','Sig. AttentionModel-MEG',...
    'RecognitionModel-MEG','Sig. RecognitionModel-MEG'},'location','northeastoutside');
axis([0 3 -0.1 0.4]);
xlabel('Time (Seconds)');
ylabel('R');
ax = get(gca,'position');
set(gca,'position',[ax(1) ax(2) 0.6 ax(4)]);
set(gcf,'position',[223         178        1329         570])