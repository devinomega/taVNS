function ECGbaselinePreProc(varargin)
%% ECG_preProc
%
% Process .cnt files and extract peaks for all blocks, pre, stim, post
%
% Input:
%       startPnt - where to start the processing (num)
%       fPath - where the files are located
%       skip - TRUE/FALSE a flag that skips over ECG files that already
%       exist
%
%
% CNT file - physiological recording from ANT
% EDF file - eyetracking data
% MAT file - basic information form the run
%
% Older files had a different file length - this accounts for this
% oldFileAdd = 120000;
% xx = (0:(135*2000-1))/2000;

%% Input Parameters
p = inputParser;

defaultPath = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\subjData\';
defaultStart= 1;
defaultSkip = true;
defaultStartName = '';

addParameter(p,'filePath',defaultPath,@ischar);
addParameter(p,'start',defaultStart,@isnumeric)
addParameter(p,'startName',defaultStartName,@ischar)
addParameter(p,'skip',defaultSkip,@islogical)

parse(p,varargin{:})

fPath = p.Results.filePath;
z= p.Results.start;
skip = p.Results.skip;
startName = p.Results.startName;

%% Parameters for analysis
%General
fFiles = dir([fPath '\TV*.cnt']);

if ~strcmp(startName,'')
    fubar = regexp({fFiles(:).name},[startName '.cnt']);
    z = find(~cellfun(@isempty, fubar));
end

%Segmenting
timeLen = 150;
samplingFreq = 2000;
dataLngth = timeLen*samplingFreq;
%AVG BPM
wndwsAvg = [60, 120];   %Time bins
wndwSz = diff(wndwsAvg,1,2);    %length of time bins
%AVG BPM 5 second bins
% wndwsAvg_5 = [[5, 15:5:120]', [15, 20:5:125]'];   %Time bins
wndwsAvg_5 = [[60:5:115]', [65:5:120]'];   %Time bins
wndwSz_5 = diff(wndwsAvg_5,1,2);    %length of time bins
%Max/Min analysis
maxDiff = 0;    %Start point for max
minDiff = 300;  %Start point for min
ppWndws = [60 120];   %Window of interst
%HRV 15 second bins
wndwsAvg_15 = [[60:15:105]', [75:15:120]'];   %Time bins
wndwSz_15 = diff(wndwsAvg_15,1,2);    %length of time bins
%% Warnings
% supress this warning about loading a variablethat's not there
warning('off', 'MATLAB:load:variableNotFound')
%%
% ==================================================================
% Process each EEG file
% ==================================================================
for i = z:numel(fFiles)
    %%
    curSubj = fFiles(i).name; %subject currently being worked on
    curMat = [fPath 'BIGMAT_' curSubj(1:end-3) 'mat'];
    
    %Check to see if ECG variable already exists for that subject
    if ~exist(curMat,'file') == 2
        rsp = questdlg(['There is no mat file for ' curSubj(1:end-4)], ...
            'Mat missing', ...
            'Next', 'Quit','Next');
        if strcmp(rsp, 'Next')
            continue
        elseif strcmp(rsp, 'Quit')
            break
        end
    else
        tmp = load(curMat,'ECG');
        if isfield(tmp, 'ECG')
            ECG=tmp.ECG;
            if isfield(tmp.ECG,'baseline')
                if skip
                    continue
                else
                    rsp = questdlg(['There''s already a baseline variable for ' curSubj(1:end-4)], ...
                        'Mat file', ...
                        'Next', 'Overwrite','Cancel','Next');
                    if strcmp(rsp, 'Next')
                        continue
                    elseif strcmp(rsp, 'Cancel') || isempty(rsp)
                        break
                    end
                end
            end
        else
            rsp = questdlg(['There is no ECG field for ' curSubj(1:end-4)], ...
                'ECG missing', ...
                'Next', 'Quit','Next');
            if strcmp(rsp, 'Next')
                continue
            elseif strcmp(rsp, 'Quit')
                break
            end
        end
    end
    
    data = pop_loadeep_v4([fPath curSubj]);   %load their eeg data
    
    if ~isempty(data.event)
        % From pre
        baseLine = [data.event(~cellfun(@isempty,regexp({data.event.type}, ...
            '111'))).latency];   %something may go wrong here if they don't have 96
        
        %Data
        curData = double(data.data(33,:)); %Raw ECG
        baseDiff = (baseLine-dataLngth);    %the number of samples missing from the start of data
        
        if baseDiff<0
            curData = [zeros(1,(abs(baseDiff)+1)) curData];
            baseLine = baseLine+abs(baseDiff)+1;
        end
        
        ECG.baseLineData = curData((baseLine-dataLngth):(baseLine-1));  %raw baseline data
        %Find the Bpm of data
        locs =peakPick(ECG.baseLineData,[curSubj(1:end-4) ' baseline']);
        if locs == 999
            break
        end
        
        if ~isnan(locs)
            secTime  = (0:(timeLen*samplingFreq)-1)/samplingFreq;
            locs = secTime(locs);   %Convert to seconds (from samples)
            ECG.baseLocsData = locs;
            % ==================================================================
            % Average BPM across entire one minute time period
            % ==================================================================
            %Pre, stim and Post segments of peak locations (including fractions
            %of beats)
            allLocs = locs >= wndwsAvg(:,1) & locs < wndwsAvg(:,2); %Find the locs (peaks) within each bin
            splLocs = allLocs; %Added beat before and after
            edgeDetect = diff(allLocs,1,2); %detect edge
            splLocs(edgeDetect == 1) = 1;  %add beat immediately before
            splLocs([false(size(splLocs,1),1), edgeDetect] == -1) = 1; % Add the beat immediatly after
            
            temp = zeros(size(wndwsAvg,1),1);
            for t = 1:size(allLocs,1)
                curLocs = locs(splLocs(t,:));
                bop = (numel(curLocs)-3) + ((curLocs(2) - wndwsAvg(t,1))/(curLocs(2)-curLocs(1))) ...
                    + ((wndwsAvg(t,2)- curLocs(end-1))/(curLocs(end)-curLocs(end-1)));        %Average BPM - partial
                temp(t) = (bop/wndwSz(t))*60;
            end
            
            ECG.baseAvgBPM= temp;
            % ==================================================================
            % Average BPM across 5 sec bins
            % ==================================================================
            allLocs_5 = locs >= wndwsAvg_5(:,1) & locs < wndwsAvg_5(:,2); %Find the locs (peaks) within each bin
            splLocs = allLocs_5; %Added beat before and after
            edgeDetect = diff(allLocs_5,1,2); %detect edge
            splLocs(edgeDetect == 1) = 1;  %add beat immediately before
            splLocs([false(size(splLocs,1),1), edgeDetect] == -1) = 1; % Add the beat immediatly after
            
            temp = zeros(size(wndwsAvg_5,1),1);
            for t = 1:size(allLocs_5,1)
                curData = locs(splLocs(t,:));
                bop = (numel(curData)-3) + ((curData(2) - wndwsAvg_5(t,1))/(curData(2)-curData(1))) ...
                    + ((wndwsAvg_5(t,2)- curData(end-1))/(curData(end)-curData(end-1)));        %Average BPM - partial
                temp(t) = (bop/wndwSz_5(t))*60;
            end
            ECG.baseAvgBPM_5 = temp;
            
            % ==================================================================
            % Median Average BPM across entire period
            % ==================================================================
            allLocs = locs >= wndwsAvg(:,1) & locs < wndwsAvg(:,2); %Find the locs (peaks) within each bin
            temp = zeros(size(wndwsAvg,1),1);
            
            for t = 1:size(allLocs,1)
                temp(t) = 60./median(diff(locs(allLocs(t,:))));
            end
            
            ECG.baseMedBPM = temp;
            
            % ==================================================================
            % Median Average BPM across 5 sec bins
            % ==================================================================
            allLocs_5 = locs >= wndwsAvg_5(:,1) & locs < wndwsAvg_5(:,2); %Find the locs (peaks) within each bin
            temp = zeros(size(wndwsAvg_5,1),1);
            
            for t = 1:size(allLocs_5,1)
                temp(t) = 60./median(diff(locs(allLocs_5(t,:))));
            end
            
            ECG.baseMedBPM_5 = temp;
            % ==================================================================
            % Max BPM for Stim
            % ==================================================================
            % Find the max BPM for the Stim time using a sliding window with 1
            % second overalp
            % time windows for stim and post (stops 10 seconds before end of
            % window
            
            curWndw = ppWndws;
            for t = curWndw(1):curWndw(2)
                %Design the window (10s) and find the locations
                window = [t t+10];
                
                Indx = find(locs >= window(1) &  locs <= window(2));
                partLoc = locs([Indx(1)-1 Indx Indx(end)+1]);  % Add the beat immediatly before and after
                part = (numel(Indx)-1) + (partLoc(2)-window(1))/(partLoc(2)-partLoc(1)) + ...
                    (window(2)-partLoc(end-1))/(partLoc(end)-partLoc(end-1));  %calculate partial beats in a period of time
                
                curBPM = (part/diff(window))*60;
                %Update current max
                if curBPM > maxDiff
                    maxTime = window;
                    maxDiff = curBPM;
                end
                if curBPM < minDiff
                    minTime = window;
                    minDiff = curBPM;
                end
            end
            ECG.baseMinBPM = [minTime minDiff];
            ECG.baseMaxBPM = [maxTime maxDiff];
            
            allLocs_15 = locs >= wndwsAvg_15(:,1) & locs < wndwsAvg_15(:,2); %Find the locs (peaks) within each bin
            % ==================================================================
            % SDRR Standard deviation of RR intervals
            % STD([RR1, RR2,...])
            % ==================================================================
            ECG.basehrVarData = std(diff(locs(allLocs(1,:))));
            
            % Fifteen seconds windows
            temp = zeros(size(wndwsAvg_15,1),1);
            for t = 1:size(allLocs_15,1)
                temp(t) = std(diff(locs(allLocs_15(t,:))));
            end
            ECG.basehrVarData_15 = temp;
            % ==================================================================
            % RMSSD - square root of the mean squared difference of
            % successive NN intervals
            % sqrt(mean([(RR1-RR2)^2,(RR2-RR3)^2...]))
            % ==================================================================
            %RMSSD
            ECG.baseRMSSD = sqrt(mean(diff(diff(locs(allLocs(1,:)))).^2));
            % Fifteen seconds windows
            temp = zeros(size(wndwsAvg_15,1),1);
            for t = 1:size(allLocs_15,1)
                temp(t) = sqrt(mean(diff(diff(locs(allLocs_15(t,:)))).^2));
            end
            ECG.baseRMSSD_15 = temp;
            % ==================================================================
            % pRR50 percent of succesive RR intervals greater than 50ms
            % ==================================================================
            ECG.basepRR50 = sum(diff(diff(locs(allLocs(1,:))))>.05)/numel(diff(diff(locs(allLocs(1,:)))));
            
            temp = zeros(size(wndwsAvg_15,1),1);
            for t = 1:size(allLocs_15,1)
                temp(t) = sum(diff(diff(locs(allLocs_15(t,:))))>.05)/numel(diff(diff(locs(allLocs(1,:)))));
            end
            ECG.basepRR50_15 = temp;
        else
            ECG.baseAvgBPM = locs;
            ECG.baseAvgBPM_5= locs;
            ECG.baseMedBPM = locs;
            ECG.baseMedBPM_5 = locs;
            ECG.baseMinBPM = locs;
            ECG.baseMaxBPM = locs;
            ECG.basehrVarData = locs;
            ECG.baseRMSSD = locs;
            ECG.basepRR50= locs;
            ECG.basehrVarData_15 = locs;
            ECG.baseRMSSD_15 = locs;
            ECG.basepRR50_15= locs;
        end
        if locs == 999
            break
        end
        %Only save the ECG data if all blocks are done - otherwise you'll have
        %partially complete data
        save(curMat, 'ECG', '-append')
    else
        uiwait(msgbox(['There is no Evts for ' curSubj(1:end-4)], 'EVT missing'));
    end
end
% Turn the warning back on
warning('off', 'MATLAB:load:variableNotFound')
