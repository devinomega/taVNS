function breathBaselinePreProc(varargin)
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

defaultPath = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData\';
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
            RespRate=tmp.RespRate;
            if isfield(tmp.RespRate,'baseRespRate')
                if skip
                    continue
                else
                    rsp = questdlg(['There''s already a baseline variable for Respiration for ' curSubj(1:end-4)], ...
                        'Mat file', ...
                        'Next', 'Overwrite','Cancel','Next');
                    if strcmp(rsp, 'Next')
                        continue
                    elseif strcmp(rsp, 'Cancel') || isempty(rsp)
                        break
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
    end
         
    
    data = pop_loadeep_v4([fPath curSubj]);   %load their eeg data
    
    if ~isempty(data.event)
        % From pre
        baseLine = [data.event(~cellfun(@isempty,regexp({data.event.type}, ...
            '111'))).latency];
        %Data
        curData = double(data.data(34,:)); %Raw RespRate
        
        blDiff = baseLine - dataLngth;
        if blDiff < 0 
            baseLine = baseLine+abs(blDiff);
        end
        RespRate.baseLineData = curData((baseLine-dataLngth):(baseLine-1));  %raw baseline data
        %% Main for loop           
            %Find the Bpm of data
            [minL, maxL] =respPick(RespRate.baseLineData,[curSubj(1:end-4) ' Baseline']);
            if minL == 999
                break
            end
            RespRate.baseLocsMin = {minL};
            RespRate.baseLocsMax= {maxL};
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
