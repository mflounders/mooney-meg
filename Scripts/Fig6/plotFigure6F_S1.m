%% commonality analysis (2 models) plotting - Fig. 6F and 6-S1

% 1) presets
% 2) loads data and relevant statistics
% 3) plots

%% plotting figure 6-figure supplement 1
% quadratic axis

ROItoPlot = [2:6 8:10 12:14 16:18 20];

%plot parameters
colors = 'rgbk';
tms = -2.6:0.01:2.9-0.01;
tytick = [0 0.01 0.02 0.03 0.04 0.06 0.08 0.1];
colors2 = 'rbk';

% close all
counter = 1;
figure
for i_ROI = ROItoPlot
    
    load(['Data/Processed/Figure6/Commonality/statistics_FMRI_MEG_Model_Commonality_' num2str(i_ROI) '.mat']);
    
    set(gcf,'position',[37         127        1416         718])
    
    subplot(5,3,counter);
    s = sign(RSquaredFmriMEG);
    area(tms,s.*sqrt(abs(RSquaredFmriMEG)),'FaceColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
        'EdgeColor',[0.501960813999176 0.501960813999176 0.501960813999176]); %,'color',colors(4));
    hold on
    pcks = find(clusterTimeSeriesRSquaredFmriMEG);
    plot(tms(pcks), 0.33 * ones(length(pcks),1),'.','color',colors(4));
    for i_model = 1:2
        s = sign(commonality2Models(i_model,:));
        plot(tms,s.*(sqrt(abs(commonality2Models(i_model,:)))),'color',colors2(i_model));
        pcks = find(clusterTimeSeriesCommonality2Models(:,i_model));
        if length(pcks) > 0
            plot(tms(pcks), (0.33-(i_model*0.03)) * ones(length(pcks),1),'.','color',colors2(i_model));
        else
            plot(0, 0,'.');
        end
    end
    legend({'MEG-FMRI', 'Sig. MEG-FMRI',...
        'Commonality Stimulus','Sig. Commonality Stimulus',...
        'Commonality Recognition','Sig. Commonality Recognition'},'location','northeastoutside');
    set(gca,'ytick',sqrt(abs(tytick)));
    axis([0 3 0 0.34]);
    set(gca,'yticklabel',tytick);
    xlabel('Time (Seconds)');
    ylabel('R^2');
    
    title(ROI_labels{i_ROI});
    
    line([0 0],[-0.05 0.34],'Color','k','LineWidth',1.5); line([2 2],[-0.05 0.34],'Color','k','LineWidth',1.5);
    counter = counter+1;
end

%% manuscript figure 5 panel F, 5 ROIs
% quadratic axis

figure; counter = 1;
for i_ROI = [1 7 11 15 19]
    
    load(['Data/Processed/Figure6/Commonality/statistics_FMRI_MEG_Model_Commonality_' num2str(i_ROI) '.mat']);
    
    tms = -2.6:0.01:2.9-0.01;
    
    set(gcf,'position',[37         127        1416         718])
    
    tytick = [0 0.01 0.02 0.03 0.04 0.06 0.08 0.1];
    
    colors2 = 'rbk';
    subplot(2,3,counter);
    s = sign(RSquaredFmriMEG);
    area(tms,s.*sqrt(abs(RSquaredFmriMEG)),'FaceColor',[0.800000011920929 0.800000011920929 0.800000011920929],...
        'EdgeColor',[0.501960813999176 0.501960813999176 0.501960813999176]);
    hold on
    pcks = find(clusterTimeSeriesRSquaredFmriMEG);
    plot(tms(pcks), 0.33 * ones(length(pcks),1),'.','color',colors(4));
    for i_model = 1:2
        s = sign(commonality2Models(i_model,:));
        plot(tms,s.*(sqrt(abs(commonality2Models(i_model,:)))),'color',colors2(i_model));
        pcks = find(clusterTimeSeriesCommonality2Models(:,i_model));
        if length(pcks) > 0
            plot(tms(pcks), (0.33-(i_model*0.03)) * ones(length(pcks),1),'.','color',colors2(i_model));
        else
            plot(0, 0,'.');
        end
    end
    legend({'MEG-FMRI', 'Sig. MEG-FMRI',...
        'Commonality Stimulus','Sig. Commonality Stimulus',...
        'Commonality Recognition','Sig. Commonality Recognition'},'location','northeastoutside');
    set(gca,'ytick',sqrt(abs(tytick)));
    axis([0 3 0 0.34]);
    set(gca,'yticklabel',tytick);
    xlabel('Time (Seconds)');
    ylabel('R^2');
    title(ROI_labels{i_ROI});
    
    counter = counter+1;
end

line([0 0],[-0.05 0.34]);
line([2 2],[-0.05 0.34]);
