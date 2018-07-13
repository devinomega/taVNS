function breath_preProc(varargin)
%% RespRate_preProc
%
% Process .cnt files and extract respiration info
%
% Input:
%       startPnt - where to start the processing (num)
%       fPath - where the files are located
%       skip - TRUE/FALSE a flag that skips over RespRate files that already
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
%Filter
sampleRate = 2000; % Hz
highEnd = 15; % Hz
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
    
    %Check to see if RespRate variable already exists for that subject
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
        tmp = load(curMat,'RespRate');
        if isfield(tmp, 'RespRate')
            if skip
                continue
            else
                rsp = questdlg(['There''s already a RespRate variable for ' curSubj(1:end-4)], ...
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

        %Data
        curData = double(data.data(34,:)); %Raw RespRate
        RespRate.locsMin = cell(1,5);   %Locations of troughs 
        RespRate.locsMax = cell(1,5);   %Locations of peaks
        RespRate.allData = cellfun(@(i) curData(i:(i+dataLngth-1)), num2cell(pre),'UniformOutput',false);
        
        % ==================================================================
        % Outcome Measures
        % ==================================================================
        %     RespRate.AvgPartBPM = zeros(5,3); % Avg BPM for [pre stim post] across a fixed window(include fractional beats)
        RespRate.min = zeros(5,6);  % Troughs of respiration rate 
        RespRate.max= zeros(5,6);   % same with min (first 3 are stim, second 3 are post)
        %% Main for loop
        
        % n is each block

        for n = 1:numel(RespRate.allData)
            % Current block of 5
            % Shorten to 135 s long
            temp = RespRate.allData{n};
             %temp = -1*RespRate.allData{n};  %in case you need to invert
            
            %Find the Bpm of data
            [minL, maxL] =respPick(temp,[curSubj(1:end-4) '_B' num2str(n)]);
            if minL == 999
                break
            end
            
            RespRate.locsMin(1,n) = {minL};
            RespRate.locsMax(1,n) = {maxL};
             findpeaks
            
        end
        if minL == 999
            break
        end
        %Only save the RespRate data if all blocks are done - otherwise you'll have
        %partially complete data
        save(curMat, 'RespRate', '-append')
    else
        uiwait(msgbox(['There is no Evts for ' curSubj(1:end-4)], 'EVT missing'));
    end
end
% Turn the warning back on
warning('off', 'MATLAB:load:variableNotFound')
