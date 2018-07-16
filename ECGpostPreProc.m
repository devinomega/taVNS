function ECGpostPreProc(varargin)
%% ECG_ECGpostPreProc
%
% Process .cnt files and extract peaks for all blocks, last post
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
timeLen = 30;
samplingFreq = 2000;
dataLngth = timeLen*samplingFreq;
%AVG BPM
wndwsAvg = [0, 30];   %Time bins
wndwSz = diff(wndwsAvg,1,2);    %length of time bins
%AVG BPM 5 second bins
wndwsAvg_5 = [[0:5:25]', [5:5:30]'];   %Time bins
wndwSz_5 = diff(wndwsAvg_5,1,2);    %length of time bins
%Max/Min analysis
maxDiff = 0;    %Start point for max
minDiff = 300;  %Start point for min
ppWndws = {0 ,30};   %Window of interst
%HRV 15 second bins
wndwsAvg_15 = [[0 15]', [15 30]'];   %Time bins
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
            if isfield(tmp.ECG,'postData')
                if skip
                    continue
                else
                    rsp = questdlg(['There''s already a post variable for ' curSubj(1:end-4)], ...
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
    
    data = pop_loadeep_v4([fPath curSubj]);   %load eeg data
    
    if ~isempty(data.event)
        % Identify post event flags
        post = [data.event(~cellfun(@isempty,regexp({data.event.type}, ...
            '1\d{1}5'))).latency];
        
        %Data
        curData = double(data.data(33,:)); %Raw ECG
        
        postDiff = numel(curData)-(post(5)+dataLngth);    %the number of samples missing from the end of data
        
        if postDiff<0
            curData = [curData zeros(1,(abs(baseDiff)+1))];
            baseLine = baseLine+abs(baseDiff)+1;
        end
        
        ECG.postData = cellfun(@(i) curData(i:(i+dataLngth-1)), num2cell(post),'UniformOutput',false);
        
         % ==================================================================
        % Outcome Measures
        % ==================================================================
        ECG.postLocsData = cell(1,5);   %Locations of peaks
        ECG.postAvgBPM = zeros(5,1); % Avg BPM for post across a fixed window
        ECG.postAvgBPM_5 = zeros(5,6); % Avg BPM for post across a fixed window
        ECG.postMedBPM = zeros(5);  % Median Average BPM across entire period
        ECG.postMedBPM_5 = zeros(5,6);    % Median Average BPM across entire period
%         ECG.postMinBPM = zeros(5,3);  % Min BPM in a sliding mindow for stim and post baseline corrected (include fractional beats)
%         ECG.postMaxBPM = zeros(5,3);   % same with min 
        ECG.posthrVarData = zeros(5); %heart rate variablity for post
        ECG.postRMSSD =  zeros(5); %Root mean square of the successive differnce post
        ECG.postpRR50 = zeros(5); %RMSSD for post
        ECG.posthrVarData_15 = zeros(5,2); %heart rate variablity for post 15 second bins
        ECG.postRMSSD_15 =  zeros(5,2); %Root mean square of the successive differnce post 15 second bins
        ECG.postpRR50_15 = zeros(5,2); %RMSSD for post 15 second bins
        %% Loop through the blocks
        for n = 1:numel(ECG.postData)
            % Current block of 5
            % Shorten to 30s long
            temp = ECG.postData{n};
            temp= temp(1:60000);
            
            %Find the Bpm of data IMPORTANT
            % This creates the window where you pick peaks from
            locs =peakPick(temp,[curSubj(1:end-4) '_B' num2str(n)]);
            %This checks for a rejected block
            if locs == 999
                break
            end
            
            %Checks for no numbers in locations
            if ~isnan(locs)
                secTime  = (0:(timeLen*samplingFreq)-1)/samplingFreq; %converted samples to seconds ona timeline
                locs = secTime(locs);   %Convert to seconds (from samples)
                ECG.postLocsData(n) = {locs};   %Make the value a structure 
                
                % ==================================================================
                % Average BPM across post
                % ==================================================================
                % Post segments of peak locations (including fractions of beats)
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
                ECG.postAvgBPM(n) = temp;
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
                ECG.postAvgBPM_5(n,:) = temp;
                
                % ==================================================================
                % Median Average BPM across entire period
                % ==================================================================
                allLocs = locs >= wndwsAvg(:,1) & locs < wndwsAvg(:,2); %Find the locs (peaks) within each bin
                temp = zeros(size(wndwsAvg,1),1);
                
                for t = 1:size(allLocs,1)
                    temp(t) = 60./median(diff(locs(allLocs(t,:))));
                end
                
                ECG.postMedBPM(n,:) = temp;
                
                % ==================================================================
                % Median Average BPM across 5 sec bins
                % ==================================================================
                allLocs_5 = locs >= wndwsAvg_5(:,1) & locs < wndwsAvg_5(:,2); %Find the locs (peaks) within each bin
                temp = zeros(size(wndwsAvg_5,1),1);
                
                for t = 1:size(allLocs_5,1)
                    temp(t) = 60./median(diff(locs(allLocs_5(t,:))));
                end
                
                ECG.postMedBPM_5(n,:) = temp;
                % ==================================================================
                % Max/Min BPM for post
                % ==================================================================
                % Find the max BPM for the Stim time using a sliding window with 1
                % second overalp
                % time windows for stim and post (stops 10 seconds before end of
                % window
%                 
%                 curWndw = ppWndws;
%                 for t = curWndw(1):curWndw(2)
%                     %Design the window (10s) and find the locations
%                     window = [t t+10];
%                     
%                     Indx = find(locs >= window(1) &  locs <= window(2));
%                     partLoc = locs([Indx(1)-1 Indx Indx(end)+1]);  % Add the beat immediatly before and after
%                     part = (numel(Indx)-1) + (partLoc(2)-window(1))/(partLoc(2)-partLoc(1)) + ...
%                         (window(2)-partLoc(end-1))/(partLoc(end)-partLoc(end-1));  %calculate partial beats in a period of time
%                     
%                     curBPM = (part/diff(window))*60;
%                     %Update current max
%                     if curBPM > maxDiff
%                         maxTime = window;
%                         maxDiff = curBPM;
%                     end
%                     if curBPM < minDiff
%                         minTime = window;
%                         minDiff = curBPM;
%                     end
%                 end
%                         ECG.postMinBPM(n,:) = [minTime minBPM];
%                     ECG.postMaxBPM(n,:) = [maxTime maxBPM];
                
                allLocs_15 = locs >= wndwsAvg_15(:,1) & locs < wndwsAvg_15(:,2); %Find the locs (peaks) within each bin
                % ==================================================================
                % SDRR Standard deviation of RR intervals
                % STD([RR1, RR2,...])
                % ==================================================================
                ECG.posthrVarData(n,:) = std(diff(locs(allLocs(1,:))));
            
                % Fifteen seconds windows
                temp = zeros(size(wndwsAvg_15,1),1);
                for t = 1:size(allLocs_15,1)
                    temp(t) = std(diff(locs(allLocs_15(t,:))));
                end
                ECG.posthrVarData_15(n,:) = temp;
                % ==================================================================
                % RMSSD - square root of the mean squared difference of
                % successive NN intervals
                % sqrt(mean([(RR1-RR2)^2,(RR2-RR3)^2...]))
                % ==================================================================
                 %RMSSD
                ECG.posteRMSSD(n,:) = sqrt(mean(diff(diff(locs(allLocs(1,:)))).^2));
                % Fifteen seconds windows
                temp = zeros(size(wndwsAvg_15,1),1);
                for t = 1:size(allLocs_15,1)
                    temp(t) = sqrt(mean(diff(diff(locs(allLocs_15(t,:)))).^2));
                end
                ECG.postRMSSD_15(n,:) = temp;
                % ==================================================================
                % pRR50 percent of succesive RR intervals greater than 50ms
                % ==================================================================
                ECG.postpRR50(n,:) = sum(diff(diff(locs(allLocs(1,:))))>.05)/numel(diff(diff(locs(allLocs(1,:)))));
            
                temp = zeros(size(wndwsAvg_15,1),1);
                for t = 1:size(allLocs_15,1)
                    temp(t) = sum(diff(diff(locs(allLocs_15(t,:))))>.05)/numel(diff(diff(locs(allLocs(1,:)))));
                end
                ECG.postpRR50_15(n,:) = temp;
            else
                ECG.postAvgBPM(n,:) = locs;
                ECG.postAvgBPM_5(n,:)= locs;
                ECG.postMedBPM(n,:) = locs;
                ECG.postMedBPM_5(n,:) = locs;
                ECG.postMinBPM(n,:) = locs;
                ECG.postMaxBPM(n,:) = locs;
                ECG.posthrVarData(n,:) = locs;
                ECG.postRMSSD(n,:) = locs;
                ECG.postpRR50(n,:) = locs;
                ECG.posthrVarData_15(n,:) = locs;
                ECG.postRMSSD_15(n,:) = locs;
                ECG.postpRR50_15(n,:) = locs;
            end
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
