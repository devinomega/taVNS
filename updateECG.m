%% Change ECG file Names to include dates
% ECG files .cnt and .evt files are saved as dates under the same subject
% We would like to save them as TV_subject#_session#
%
% List all the .mat files
% Grab their date
% Look for that date with that subject
% Rename said file
%If you can't find the corresponding file - a dialogue box will pop up

%%Path of Subject data
fpath = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData';
sfiles = dir([fpath '\TV*.mat']);

efiles = dir([fpath '\TV*.cnt']);
eNames = {efiles(:).name};

%%
for i = 1:numel(sfiles)
    %extract the session and subject ID from the *.mat name
    curSubj = sfiles(i).name;
    session = regexp(curSubj,'TV_(\d\d)_(\d\d).mat','tokens','once');
    
    %Check to make sure there is an event file
    if ~(exist([fpath '\' curSubj(1:end-3) 'evt'],'file') == 2)
         % load run_date
        load([fpath '\' curSubj],'run_Date'); %load data
        ecgDate = datestr(run_Date,'yyyy-mm-dd');
        
        %Check to make sure the data in the ECG and the run data are the
        %same
        a = ~cellfun(@isempty,regexp(eNames,['[a-zA-Z]*_' session{1} '(_\d\d)?_[a-zA-Z]*_' ecgDate '_*']));
        if sum(a) ~= 1
           answer = questdlg(['Can''t finding a matching date. Or file' ...
               ' does not exist for ' curSubj(1:end-4)], 'Date error', ...
                'OK, continue','Break','OK, continue');
            if strcmp(answer, 'Break')
                return
            end
        else
            cntF = [fpath '\' eNames{a}];
            evtF = [fpath '\' eNames{a}(1:end-3) 'evt'];
            
            %Now find the matching ECG file and date
            if ~(exist(cntF,'file') == 2) || ~(exist(evtF,'file') == 2)
                answer = questdlg(['Can''t finding a matching ECG File'...
                    curSubj(1:end-4)], 'ECG error', ...
                    'OK, continue','Break','OK, continue');
                if strcmp(answer, 'Break')
                    return
                end
            else
                % Otherwise - change the names to the "TV_(\d\d)_(\d\d)"
                % format
                movefile(cntF,[fpath '\' curSubj(1:end-3) 'cnt']);
                movefile(evtF,[fpath '\' curSubj(1:end-3) 'evt']);
            end

        end
        
    end
end

