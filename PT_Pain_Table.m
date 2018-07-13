fPath = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData\';  %subject Files
fFiles = dir([fPath '\BIGMAT*.mat']);
numSubj = regexp({fFiles(:).name},'BIGMAT_TV_0(\d)_*|BIGMAT_TV_([1-9]\d)_*','tokens','once');
numSubj = max(unique(str2double([numSubj{:}])));
analyses = {'PT', 'Dose','Pain','Pre'};
ordr = {'Active_25 Hz','Control_25 Hz', 'Active_100 Hz', 'Control_100 Hz', ...
    'Active_Burst', 'Control_Burst', 'Active_Bilateral 25 Hz', 'Control_Bilateral 25 Hz'};

dataMat = NaN(numel(ordr),numel(analyses),numSubj);
for i = 1:numel(fFiles)
    curSubj = fFiles(i).name;
    curData = load([fPath curSubj],'sessNum','pain','stimInt');
    
    subjNum = regexp(curSubj,'BIGMAT_TV_(\d\d)_*','tokens','once');
    subjNum = str2double(subjNum{:});
    curSess = strcmp(curData.sessNum,ordr); %Identify the session 
    
    dataMat(curSess,1:2,subjNum) = curData.stimInt(1,[1,3]);

    dataMat(curSess,3,subjNum) = nanmean(curData.pain(2:6));
    dataMat(curSess,4,subjNum) = curData.pain(1);
end
dataMat(dataMat == 999) = NaN;

nanmean(dataMat,3)

x = squeeze(dataMat(7,4,:));
y = squeeze(dataMat(8,4,:));
[h, p] = ttest(x,y);