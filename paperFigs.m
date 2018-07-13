%% Figures/ANOVA 

fFiles = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData\';  %subject Files
ordr = {'Active_25 Hz','Control_25 Hz', 'Active_100 Hz', 'Control_100 Hz', ...
    'Active_Burst', 'Control_Burst', 'Active_Bilateral 25 Hz', 'Control_Bilateral 25 Hz'};
fiveSec = {'pre','stim_05','stim_10','stim_15','stim_20','stim_25','stim_30','stim_35','stim_40','stim_45','stim_50','stim_55','stim_60'};
rmv = zeros([29 1]);
rmv([1,7,16,21]) = 1;
%% Plot Active vs Sham, collapsed across block, use 5 second bins, not baseline corrected
curData = load([fFiles 'medBPM_5_Group_Block.mat']);
curData = curData.strctData;
curData(logical(rmv),:,:,:) = [];

groupings = [[1,3,5,7];[2,4,6,8]];

figure
    hold on
for f= 1:size(groupings,1)
    slctData = reshape(permute(curData(:,:,groupings(f,:)),[2,3,1]),[27,100]);
    SEM = std(slctData,0,2,'omitnan')/sqrt(size(slctData,2));
    plotData = nanmean(nanmean(nanmean(curData(:,:,:,groupings(f,:)),4),3),1);
%     errorbar(plotData',SEM','-o','Linewidth',4)
    plot(plotData,'-o','Linewidth',4);

%     plot([0 13], [0 0],'--k')
end

box off

ttl = 'Active vs. Sham'; 
title(ttl, 'FontSize', 20)
xticks(1:size(plotData,2))
ylabel('BPM')
xlabel('Time (sec)')
% xlim([.75 12.25])
% ylim([-2 3])
xticklabels(5:5:135)
% hax=axes;  
line([3 3],ylim,'Color','black','LineStyle','--')
line([15 15],ylim,'Color','black','LineStyle','--')
legend({'Active', 'Sham','Start Stim','End Stim'},'Location','Northwest')
hold off

%% Plot each dose (active vs sham), collapsed across block, use 5 second bins, not baseline corrected
curData = load([fFiles 'medBPM_5_Group_Block.mat']);
curData = curData.strctData;
curData(logical(rmv),:,:,:) = [];
groupings = [1:2;3:4;5:6;7:8];

figure
hold on
for m =1:size(groupings,1)
    subplot(2,2,m)
    ttl = regexp(ordr(groupings(m,1)),'.*_(.*)','tokens','once');
    ttl = ttl{:};
    slctData = reshape(permute(curData(:,:,:,groupings(m,:)),[4,2,3,1]),[2,27,125]);
    SEM = std(slctData,0,3,'omitnan')/sqrt(size(slctData,3));     
    plotData = squeeze(nanmean(nanmean(curData(:,:,:,groupings(m,:)),3),1));
    plot(plotData,'-o','Linewidth',4);
%     errorbar(plotData,SEM','-o','Linewidth',4)
    title(ttl{:}, 'FontSize', 20)
    xticks(1:size(plotData,1))
    xticklabels(5:5:135)
    ylabel('BPM')
    xlabel('Time (sec)')
    box off
    line([3 3],ylim,'Color','black','LineStyle','--')   
    line([15 15],ylim,'Color','black','LineStyle','--')
    legend([strrep(ordr(groupings(m,:)),'_', ' '), 'Start Stim', 'End Stim']','Location','Southeast')
end

%% Plot each dose (active vs sham), collapsed across block, use 5 second bins, not baseline corrected
curData = load([fFiles 'medBPM_5_Group_Block_BLC.mat']);
curData = curData.strctData;
curData(logical(rmv),:,:,:) = [];
groupings = [1:2;3:4;5:6;7:8];

figure
hold on
for m =1:size(groupings,1)
    subplot(2,2,m)
    ttl = regexp(ordr(groupings(m,1)),'.*_(.*)','tokens','once');
    ttl = ttl{:};
    slctData = reshape(permute(curData(:,:,:,groupings(m,:)),[4,2,3,1]),[2,size(curData,2),125]);
    SEM = std(slctData,0,3,'omitnan')/sqrt(size(slctData,3));     
    plotData = squeeze(nanmean(nanmean(curData(:,:,:,groupings(m,:)),3),1));
    plot(plotData,'-o','Linewidth',4);
%     errorbar(plotData,SEM','-o','Linewidth',4)
    hold on
    title(ttl{:}, 'FontSize', 20)
    xticks(1:size(plotData,1))
    xticklabels(15:5:135)
    ylabel('Change in BPM from Baseline')
    xlabel('Time (sec)')
    box off
    line(xlim, [0 0],'Color','black','LineStyle','--')  
    line([15 15],ylim,'Color','black','LineStyle','--')
    legend([strrep(ordr(groupings(m,:)),'_', ' '),  'End Stim']','Location','Southeast')
end


%% all five blocks 
% Not doing the SEM correct (different for each) but close enough
curData = load([fFiles 'medBPM_5_Group_Block_BLC.mat']);
curData = curData.strctData;
curData(logical(rmv),:,:,:) = [];
groupings = [1:2;3:4;5:6;7:8];

figure
hold on
for m =1:numel(ordr)
    subplot(4,2,m)
    ttl = strrep(ordr{m},'_', ' ');
    temp = curData(:,:,:,m);
    
    SEM = squeeze(std(temp,0,1,'omitnan'))/sqrt(size(temp,1));     
    plotData = squeeze(nanmean(temp,1));
    plot(plotData,'-o','Linewidth',4);
%     errorbar(plotData,SEM,'-o','Linewidth',4)
    hold on
    title(ttl, 'FontSize', 20)
    xticks(1:size(plotData,1))
    xticklabels(15:5:135)
    ylabel('Change in BPM from Baseline')
    xlabel('Time (sec)')
    box off
    line(xlim, [0 0],'Color','black','LineStyle','--')  
    line([15 15],ylim,'Color','black','LineStyle','--')
    legend({'Block 1','Block 2','Block 3','Block 4','Block 5','End Stim'},'Location','Southeast')
    xlim([0 27])
end

%% HRV 15 second bins 
curData = load([fFiles 'RMSSD_15_Group.mat']);
% curData = load([fFiles 'hrVarData_Group_BLC.mat']);
curData = curData.strctData;
curData(logical(rmv),:,:,:) = [];
actv = curData(:,:,[1:2:8]);
sham = curData(:,:,[2:2:8]);
ttl = {'Pre','Stim','Recovery'};

SEM = zeros(9,2);

for n = 1:9
    SEM(n,1) = std(reshape(squeeze(actv(:,n,:)), [100 1]),'omitnan')/sqrt(numel(reshape(squeeze(actv(:,n,:)), [100 1]))) ;
    SEM(n,2) = std(reshape(squeeze(sham(:,n,:)), [100 1]),'omitnan')/sqrt(numel(reshape(squeeze(sham(:,n,:)), [100 1]))); 
end
AVG = zeros(9,2);
AVG(:,1) = nanmean(nanmean(actv,1),3)';
AVG(:,2) = nanmean(nanmean(sham,1),3)';

bns = [[1,2,6];[1,5,9]];

figure

for b = 1:size(bns,2)
    subplot(1,3,b)
    pltData = AVG(bns(1,b):bns(2,b),:);
%     errBrs = SEM(bns(1,b):bns(2,b),:);
    bar(pltData,'BarWidth', 1)
%         hold on
%     errorbar(pltData,errBrs,'.')
    title(ttl{b}, 'FontSize', 20)
    legend show
    legend({'Active','Sham'})
end

% actv = permute(curData(:,:,[1:2:8]),[2 1 3]);

% SEMA = std(reshape(actv,[25 36]),0,3,'omitnan')/sqrt(size(rshp,2));
% rshp = reshape(permute(curData,[3 1 2]),[8 25*9]);
% SEM= std(curData,0,3,'omitnan')/sqrt(size(rshp,2));
avgACT= squeeze(nanmean(nanmean(actv,3),1));
avgSham= squeeze(nanmean(nanmean(sham,3),1));
plotData = [avgACT;avgSham];

figure
bar(plotData(:,1),'BarWidth', 1)
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

