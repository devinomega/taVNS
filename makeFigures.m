%% Figures/ANOVA 

fFiles = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData\';  %subject Files
ordr = {'Active_25 Hz','Control_25 Hz', 'Active_100 Hz', 'Control_100 Hz', ...
    'Active_Burst', 'Control_Burst', 'Active_Bilateral 25 Hz', 'Control_Bilateral 25 Hz'};
fiveSec = {'pre','stim_05','stim_10','stim_15','stim_20','stim_25','stim_30','stim_35','stim_40','stim_45','stim_50','stim_55','stim_60'};

% 
% t = array2table([[1:25]',dataShape],'VariableNames', ['Subj',meas']);
% Meas = array2table(wiMat,'VariableNames',{'Condition', 'Block','PreStim'});
% rm = fitrm(t,'Active_25_Hz_01_pre-Control_Bilateral_25_Hz_05_stim~Subj','WithinDesign',Meas);
% 
% 
% grp = [61 80];
% t = array2table([[1:25]',dataShape(:,grp)],'VariableNames', ['Subj',meas(grp)']);
% Meas = array2table(wiMat(grp,:),'VariableNames',{'Condition', 'Block','PreStim'});
% rm = fitrm(t,[meas{grp(1)} '-' meas{grp(2)} '~Subj'],'WithinDesign',Meas);
% 
% [ranovatbl,A,C,D] = ranova(rm,'WithinModel','Condition+Block+PreStim');
% 
% 
% % figure
% for m =1:8
%     for n = 1:5
%         subplot(8,5,n+(m-1)*5)
%         plot(permute(squeeze(allData(m,:,n,:)),[2 1]),'-o')
%         xticks([1 2])
%         xlim([.75 2.25])
%         xticklabels({'pre' 'post'})
%     end
% end


%% Plot all data as scatter/box plot
% flag = 1; % any datapoint outside of 10 is plotted as 10
% curData = load([fFiles 'AvgBPM_5_Group_BLC.mat']);
% curMeas = strrep(ordr,'_',' '); 
% figure
% dataOrder = [1,3,5,7,2,4,6,8];
% 
% for m =1:size(curData.strctData,4)
%     subplot(2,4,m)
%     dataSet = squeeze(curData.strctData(:,:,:,dataOrder(m)));
%     dataSet = reshape(permute(dataSet,[2 1 3]),[12,29*5])';
%     
%     if flag
%         % for viewability - set a ceiling and floor
%         dataSet(dataSet>10) = 10;
%         dataSet(dataSet<-10) = -10;
%     end
%     
%     title(curMeas(dataOrder(m)))
%     hold on
%     scalerV = 5;
%     
%     for p = 1:size(dataSet,2)
%         x=ones(length(dataSet),1).*(1+(rand(length(dataSet),1)-0.5)/scalerV);
%         f1=scatter(x.*p,dataSet(:,p),'r','filled');f1.MarkerFaceAlpha = 0.4;hold on
%         scalerV = scalerV+5;
%     end
%     
%     boxplot(dataSet,'Notch','on','Labels', {'5','10','15','20','25','30','35','40','45','50','55','60'},'Whisker',1)
%     lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
%     set(lines, 'Color', 'k','Linewidth',2);
%     axesT = findobj(gcf, 'type', 'axes');
% 	set(axesT,'YLim',[-12,12])
% end
% boop =permute(curData,[2 3 4 1]);
% bigNum = find(boop > 20 | boop < -20);
% 
% subjMat = reshape(repmat(1:29,5*8*12,1),8*5*12*29,1);
% ordrMat = reshape(repmat(1:8,5*12,29),8*5*12*29,1);    %group within subje t matrix
% blockMat = reshape(repmat(1:5,12,8*29),8*5*12*29,1); %block headers
% segMat = reshape(repmat(1:12,1,5*8*29),8*5*12*29,1);
% 
% subjMat(bigNum)
% ordrMat(bigNum)
%% Plot averaged BPM data 
curData = load([fFiles 'medBPM_5_Group_BLC.mat']);
curData = curData.strctData;

groupings = [1:2;3:4;5:6;7:8];

figure
for m =1:size(groupings,1)
    subplot(2,2,m)
    ttl = regexp(ordr(groupings(m,1)),'.*_(.*)','tokens','once');
    ttl = ttl{:};
    slctData = reshape(permute(curData(:,:,:,groupings(m,:)),[4,2,3,1]),[2,12,145]);
    
    
    SEM = std(slctData,0,3,'omitnan')/sqrt(size(slctData,3));     
    plotData = squeeze(nanmean(squeeze(nanmean(curData(:,:,:,groupings(m,:)),3)),1))';
    
    errorbar(plotData',SEM','-o','Linewidth',4)
    hold on
    plot([0 13], [0 0],'--k')
    box off
    legend(strrep(ordr(groupings(m,:)),'_', ' ')','Location','Northwest')
    title(ttl{:}, 'FontSize', 20)
%     title([ttl{:} ' Stimulation Delta HR'], 'FontSize', 20)
    xticks(1:size(plotData,2))
    xlim([.75 12.25])
    ylim([-2 3])
    xticklabels(5:5:60)
    hold off
end

%%
curData = load([fFiles 'AvgBPM_5_Group.mat']);
curData = curData.strctData;

groupings = [1:2;3:4;5:6;7:8];

figure
for m =1:size(groupings,1)
    subplot(2,2,m)
    ttl = regexp(ordr(groupings(m,1)),'.*_(.*)','tokens','once');
    ttl = ttl{:};
    slctData = reshape(permute(curData(:,:,:,groupings(m,:)),[4,2,3,1]),[2,12,145]);
    
    
    SEM = std(slctData,0,3,'omitnan')/sqrt(size(slctData,3));     
    plotData = squeeze(nanmean(squeeze(nanmean(curData(:,:,:,groupings(m,:)),3)),1))';
    
    errorbar(plotData',SEM','-o','Linewidth',4)
    hold on
    plot([0 13], [0 0],'--k')
    box off
    legend(strrep(ordr(groupings(m,:)),'_', ' ')','Location','Northwest')
    title(ttl{:}, 'FontSize', 20)
%     title([ttl{:} ' Stimulation Delta HR'], 'FontSize', 20)
    xticks(1:size(plotData,2))
    xlim([.75 12.25])
    ylim([-2 3])
    xticklabels(5:5:60)
    hold off
end


%% Plot HR Variability 

curData = load([fFiles 'RMSSD_Group_BLC.mat']);
curData = load([fFiles 'hrVarData_Group_BLC.mat']);
curData = squeeze(curData.strctData);

rshp = reshape(permute(curData,[3 1 2]),[8 29*5]);
% std(rsp,0,2,'omitnan')/sqrt(size(rshp,2));
avgData = squeeze(nanmean(nanmean(curData,2),1));
plotData = reshape(avgData,[2 4])';

figure
bar(plotData,'BarWidth', 1)
xticklabels({'25 Hz','100 Hz', 'Burst', 'Bilateral'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
legend show
legend({'Active','Sham'})
for m =1:size(groupings,1)
    subplot(2,2,m)
    ttl = regexp(ordr(groupings(m,1)),'.*_(.*)','tokens','once');
    ttl = ttl{:};
    slctData = reshape(permute(curData(:,:,:,groupings(m,:)),[4,2,3,1]),[2,12,145]);

    hold on
    plot([0 13], [0 0],'--k')
    legend(strrep(ordr(groupings(m,:)),'_', ' ')','Location','Northwest')
    title(ttl{:}, 'FontSize', 20)
    box off
%     title([ttl{:} ' Stimulation Delta HR'], 'FontSize', 20)
    xticks(1:size(plotData,2))
    xlim([.75 12.25])
    ylim([-1.75 2])
    xticklabels(5:5:60)
    hold off
end

%% Plot each block around the average 
curData = load([fFiles 'hrVarData_Group_BLC.mat']);
curData = curData.strctData;
lncolor = get(gca,'colororder');
close gcf

groupings = [1:2;3:4;5:6;7:8];

figure
for m =1:size(groupings,1)
    subplot(2,2,m)
    ttl = regexp(ordr(groupings(m,1)),'.*_(.*)','tokens','once');
    ttl = ttl{:};
    
    hold on

    for g = 1:size(groupings,2)
          tWeight = 1;    %
        for b = 1:size(curData,3)
            p1 = plot(squeeze(nanmean(curData(:,:,b,groupings(m,g)),1)),'Linewidth', 3);
            set(p1,{'color'},num2cell(lncolor(g,:),2))
            p1.Color(4) = tWeight;

            tWeight = tWeight*.75;
        end
    end
    
    slctData = reshape(permute(curData(:,:,:,groupings(m,:)),[4,2,3,1]),[2,12,145]);
    SEM = std(slctData,0,3,'omitnan')/sqrt(size(slctData,3));     
    plotData = squeeze(nanmean(squeeze(nanmean(curData(:,:,:,groupings(m,:)),3)),1))';
    
    errorbar(plotData',SEM','-o')
    hold on
    plot([0 13], [0 0],'--k')
    legend(strrep(ordr(groupings(m,:)),'_', ' ')','Location','Northwest')
    title([ttl{:} ' Stimulation Delta HR'])
    xticks(1:size(plotData,2))
    xlim([.75 12.25])
    ylim([-1.75 2])
    xticklabels(5:5:60)
    hold off
end


% make a figure
% grpTime = squeeze(nanmean(subjGrpTime,2));
% figure
% title('Bilateral 25 Hz BPM Baseline Corrected')
% plot(grpTime(7:8,:)','-o')
% hold on
% plot(1:12,zeros(1,12),'--k')
% legend('Active', 'Sham','Baseline','Location','Northwest')
% xlabel('time points (5 sec bins)')
% xlim([1 12])
% xticks([5:5:60])
% ylabel('BPM change from baseline')
