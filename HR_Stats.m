function HR_Stats(varargin)
% Stats Analysis for HR tVANS
%
%   'analysisType' = cell of strings where each string is an type of
%           analysis to be done, e.g. 
%           'AVG_BPM'
%           'AvgBPM_5'
%           'hrVarData'
%           'RMSSD'
%           'PRR50'
%   'Block' - a logic array determing if you want to include or collapse
%           across blocks
%   'PreStim' - logic array determining if you want to collpase across pre
%           and the other bins
%% Input Parameters
p = inputParser;

defaultFPath = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData\';  %subject Files
defaultSPath = 'C:\Users\devinomega\Dropbox\GC\Marom\CranialNerves\Vagus\Analysis\';    %save path
defaultanalysisType= {'AvgBPM','AvgBPM_5','medBPM_5','hrVarData','RMSSD','pRR50','hrVarData_15','RMSSD_15','pRR50_15'};   % default do all analysis
defaultBlock = logical(ones([1,numel(defaultanalysisType)]));  %collapse blocks = 0; keep blocks = 1
defaultpreStim= logical(zeros([1,numel(defaultanalysisType)]));  %Keep (0) or collapse (1) baseline(5 sec bins) 

addParameter(p,'filePath',defaultFPath,@ischar);
addParameter(p,'savePath',defaultSPath,@ischar);
addParameter(p,'analysisType',defaultanalysisType,@iscell)
addParameter(p,'Block',defaultBlock,@islogical)
addParameter(p,'preStim',defaultpreStim,@islogical)

parse(p,varargin{:})

fPath = p.Results.filePath;
sPath = p.Results.savePath;
analyses = p.Results.analysisType;
incBlock = p.Results.Block;
incPre= p.Results.preStim;

%
fFiles = dir([fPath '\BIGMAT*.mat']);
numSubj = regexp({fFiles(:).name},'BIGMAT_TV_0(\d)_*|BIGMAT_TV_([1-9]\d)_*','tokens','once');
numSubj = max(unique(str2double([numSubj{:}])));
ordr = {'Active_25 Hz','Control_25 Hz', 'Active_100 Hz', 'Control_100 Hz', ...
    'Active_Burst', 'Control_Burst', 'Active_Bilateral 25 Hz', 'Control_Bilateral 25 Hz'};
blck = {'1','2','3','4','5'};
%There's a better way to do this but I am too tired
fiveSec = {'pre_05','pre_10', 'pre_15','stim_05','stim_10','stim_15','stim_20','stim_25','stim_30','stim_35','stim_40','stim_45','stim_50','stim_55','stim_60','post_05','post_10','post_15','post_20','post_25','post_30','post_35','post_40','post_45','post_50','post_55','post_60'};
fifteenSec = {'pre_15','stim_15','stim_30','stim_45','stim_60','post_15','post_30','post_45','post_60'};
oaMean = {'pre','stim','post'};


% Set up data data structures
matchAnalyses = regexp(analyses,'.*_(\d)*','tokens','once');    % How long the data is 
% mLog = ~cellfun(@isempty, matchAnalyses);

notFive = cellfun(@isempty,matchAnalyses);
prePost = cell(1,numel(analyses));
allData = cell(1,numel(analyses));
for a = 1:numel(matchAnalyses)
    if isempty(matchAnalyses{a})
        prePost(a) = {oaMean};    
        allData(a) = {nan(numSubj,numel(oaMean),numel(blck),numel(ordr))};
    elseif str2num(matchAnalyses{a}{1})==5
        prePost(a) = {fiveSec};
        allData(a) = {nan(numSubj,numel(fiveSec),numel(blck),numel(ordr))};
    elseif str2num(matchAnalyses{a}{1})==15
        prePost(a) = {fifteenSec};
        allData(a) = {nan(numSubj,numel(fifteenSec),numel(blck),numel(ordr))};
    end  

end
%% Compile Data
for i = 1:numel(fFiles)
    curSubj = fFiles(i).name;
    curData = load([fPath curSubj],'sessNum','ECG');
    
    if isfield(curData, 'ECG')
        subjNum = regexp(curSubj,'BIGMAT_TV_(\d\d)_*','tokens','once');
        curSess = strcmp(curData.sessNum,ordr); %Identify the session the data belongs to
        
        for a = 1:numel(analyses)
            curAnalyses = allData{a}; %(subjId, pre/post,block, session)
            curAnalyses(str2double(subjNum{:}),:,:,curSess) = curData.ECG.(analyses{a})(:,1:size(curAnalyses,2))'; %FLAG
            allData{a} = curAnalyses;
        end
    end
end
%% Complete Segment

% Make a table of actual data
% Make a table of the measures
% ==================================================================
% GroupXBlockXPrePost
% ==================================================================

for r = 1:numel(analyses)
    tmpSeg = prePost{r};    %segment headers
    BLC_ttl = '';   %baseline corrected
    block_ttl = ''; %block title
    
    curAnalyses = allData{r};
    %Check to see if the pre condition will be subtracted
    if ~incPre(r)    
        if str2num(matchAnalyses{a}{1})==5
            tmpSeg(1:3) = []; %remove pre from headers
            curAnalyses = curAnalyses(:,4:end,:,:)-curAnalyses(:,3,:,:);%Change data to reflect removing pre FLAG
        else
            tmpSeg(1) = []; %remove pre from headers
            curAnalyses = curAnalyses(:,2:end,:,:)-curAnalyses(:,1,:,:);%Change data to reflect removing pre FLAG
        end
        tmpSeg = strcat('BLC_', tmpSeg);    %Add the new baseline corrected tag
        BLC_ttl = '_BLC';
    end
    
    %Determine the number of total condtions
    curOrdr = size(curAnalyses,4);    %Number of conditions FLAG
    curBlock = ((size(curAnalyses,3)-1)*incBlock(r))+1;    %How many blocks
    curSeg = size(curAnalyses,2); %How many 'bins'  FLAG
    
    %Form headers
    ordrRep = reshape(repmat(ordr,curBlock*curSeg,1),curOrdr*curBlock*curSeg,1);    %group headers
    ordrMat = reshape(repmat(1:curOrdr,curBlock*curSeg,1),curOrdr*curBlock*curSeg,1);    %group within subje t matrix
    curMeas = ordrRep;
    curWiMat = ordrMat;
    
    %Determine if blocks are included
    if incBlock(r) 
        blockRep = reshape(repmat(strcat('0',blck),curSeg,curOrdr),curOrdr*curBlock*curSeg,1); %block headers
        blockMat = reshape(repmat(1:curBlock,curSeg,curOrdr),curOrdr*curBlock*curSeg,1); %block headers
        curMeas = [curMeas blockRep];
        curWiMat = [curWiMat blockMat];
        block_ttl = '_Block';
    else
        curAnalyses = squeeze(nanmean(curAnalyses,3));
    end
    
    segRep = repmat(tmpSeg',curBlock*curOrdr,1);
    segMat = repmat([1:curSeg]',curBlock*curOrdr,1);
    curWiMat = [curWiMat segMat];
    curMeas = [curMeas segRep];
    
    curMeas = join(curMeas);
    curMeas = strrep(curMeas,' ','_'); 
    
    %Reshape the data to be 2 dimensional
    %remove any subject that is missing all lines of data 
    %Add 999 as a misisng number
    dataShape = reshape(curAnalyses,[numSubj numel(curMeas)]);
    dataShape = dataShape(~all(isnan(dataShape),2),:);
    dataShape(isnan(dataShape)) = 999;

    ttl = [analyses{r} '_Group' block_ttl BLC_ttl ];
    sData = [curMeas';num2cell(dataShape)]; %data 'flattened'
    strctData = curAnalyses; %Data with a structure
    
    xlswrite([sPath ttl '.xlsx'],sData);
    save([fPath ttl '.mat'],'sData','strctData')
end


