function ECG_preProc(varargin)
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
dataLngth = 135*2000;
%AVG BPM
wndwsAvg = [[5;15;75], [15;75;135]];   %Time bins
wndwSz = diff(wndwsAvg,1,2);    %length of time bins
%AVG BPM 5 second bins
wndwsAvg_5 = [[0:5:130]', [5:5:135]'];   %Time bins
wndwSz_5 = diff(wndwsAvg_5,1,2);    %length of time bins
%Max/Min analysis
maxDiff = 0;    %Start point for max
minDiff = 300;  %Start point for min
ppWndws = {[15 75] [75 115]};   %Window of interst
%HRV 15 second bins
wndwsAvg_15 = [[0:15:130]', [15:15:135]'];   %Time bins
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
            if skip
                continue
            else
                rsp = questdlg(['There''s already a ECG variable for ' curSubj(1:end-4)], ...
                    'Mat file', ...
                    'Next', 'Overwrite','Cancel','Next');
                if strcmp(rsp, 'Next')
                    continue
                elseif strcmp(rsp, 'Cancel') || isempty(rsp)
                    break
                end
            end
        end
    end
    
    data = pop_loadeep_v4([fPath curSubj]);   %load their eeg data
    
    if ~isempty(data.event)
        % From pre
        pre = [data.event(~cellfun(@isempty,regexp({data.event.type}, ...
            '1\d{1}1'))).latency];
        
        %     Adds 60 seconds to the end of the data
        %         if isempty(postEnd)
        %             postEnd = post + 2000*60;
        %         end
        
        %Data
        curData = double(data.data(33,:)); %Raw ECG
        ECG.locsData = cell(1,5);   %Locations of peaks
        ECG.allData = cellfun(@(i) curData(i:(i+dataLngth-1)), num2cell(pre),'UniformOutput',false);
        
        % ==================================================================
        % Outcome Measures
        % ==================================================================
        %     ECG.AvgPartBPM = zeros(5,3); % Avg BPM for [pre stim post] across a fixed window(include fractional beats)
        ECG.AvgBPM = zeros(5,3); % Avg BPM for [pre stim post] across a fixed window
        ECG.AvgBPM_5 = zeros(5,27); % Avg BPM for [pre stim post] across a fixed window
        ECG.minBPM = zeros(5,6);  % Min BPM in a sliding mindow for stim and post baseline corrected (include fractional beats)
        ECG.maxBPM = zeros(5,6);   % same with min (first 3 are stim, second 3 are post)
        ECG.hrVarData = zeros(5,3); %heart rate variablity for [pre stim post]
        ECG.RMSSD =  zeros(5,3); %Root mean square of the successive differnce [pre stim post]
        ECG.pRR50 = zeros(5,3); %RMSSD for [pre stim post]
        %% Main for loop
        
        % n is each block
        for n = 1:numel(ECG.allData)
            % Current block of 5
            % Shorten to 135 s long
            temp = ECG.allData{n};
            
            %In case you need to  invert...
%             temp = temp*-1;
            
            %Find the Bpm of data
            locs =peakPick(temp,[curSubj(1:end-4) '_B' num2str(n)]);           
            if locs == 999
                break
            end
            
            if ~isnan(locs)
                secTime  = (0:(135*2000)-1)/2000;
                locs = secTime(locs);   %Convert to seconds (from samples)
                ECG.locsData(n) = {locs};
                % ==================================================================
                % Average BPM across all periods (pre, stim,post)
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
                
                ECG.AvgBPM(n,:) = temp;
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
                ECG.AvgBPM_5(n,:) = temp;
                % ==================================================================
                % Max/Min BPM for Stim
                % ==================================================================
                % Find the max BPM for the Stim time using a sliding window with 1
                % second overalp
                % time windows for stim and post (stops 10 seconds before end of
                % window
                for w = 1:numel(ppWndws)
                    curWndw = ppWndws{w};
                    for t = curWndw(1):curWndw(2)
                        %Design the window (10s) and find the locations
                        window = [t t+10];
                        
                        Indx = find(locs >= window(1) &  locs <= window(2));
                        partLoc = locs([Indx(1)-1 Indx Indx(end)+1]);  % Add the beat immediatly before and after
                        prePart = (numel(Indx)-1) + (window(1)-partLoc(2))/(partLoc(2)-partLoc(1)) + ...
                            (window(2)-partLoc(end-1))/(partLoc(end)-partLoc(end-1));  %calculate partial beats in a period of time
                        
                        curBPM = (prePart/10)*60;
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
                    ind = (w-1)*3+1:(w-1)*3+3;
                    minBPM = minDiff - ECG.AvgBPM(n,1) ;  %baseline correct
                    maxBPM = maxDiff - ECG.AvgBPM(n,1);  %baseline correct
                    ECG.minBPM(n,ind) = [minTime minBPM];
                    ECG.maxBPM(n,ind) = [maxTime maxBPM];
                end
                
                % ==================================================================
                % SDRR Standard deviation of RR intervals
                % STD([RR1, RR2,...])
                % ==================================================================
                %                 ECG.hrVarData(n,:) = [std(1./diff(locs(preLog))) std(1./diff(locs(stimLog)))...
                %                     std(1./diff(locs(postLog)))];
                
                ECG.hrVarData(n,:) = [std(diff(locs(allLocs(1,:)))) std(diff(locs(allLocs(2,:))))...
                    std(diff(locs(allLocs(3,:))))];
                % ==================================================================
                % RMSSD - square root of the mean squared difference of
                % successive NN intervals
                % sqrt(mean([(RR1-RR2)^2,(RR2-RR3)^2...]))
                % ==================================================================
                %RMSSD
                ECG.RMSSD(n,:) = [sqrt(mean(diff(diff(locs(allLocs(1,:)))).^2)) ...
                    sqrt(mean(diff(diff(locs(allLocs(2,:)))).^2)) ...
                    sqrt(mean(diff(diff(locs(allLocs(3,:)))).^2))];
                
                % ==================================================================
                % pRR50 percent of succesive RR intervals greater than 50ms
                % ==================================================================
                ECG.pRR50(n,:) = [sum(diff(diff(locs(allLocs(1,:))))>.05)/numel(diff(diff(locs(allLocs(1,:))))) ...
                    sum(diff(diff(locs(allLocs(2,:))))>.05)/numel(diff(diff(locs(allLocs(2,:))))) ...
                    sum(diff(diff(locs(allLocs(3,:))))>.05)/numel(diff(diff(locs(allLocs(3,:)))))];
            else
                ECG.AvgBPM(n,:) = locs;
                ECG.AvgBPM_5(n,:) = locs;
                ECG.minBPM(n,:) = locs;
                ECG.maxBPM(n,:) = locs;
                ECG.hrVarData(n,:) = locs;
                ECG.RMSSD(n,:)= locs;
                ECG.pRR50(n,:)= locs;
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
