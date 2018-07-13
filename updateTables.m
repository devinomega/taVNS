%create an OLE automation server for excel
e = actxserver('Excel.Application');

%% Path for demographics file
%sheet1 = Demographics and PT; sheet2 = NAS; sheet3 = randomization
% fpath = 'C:\Users\dadair\Dropbox\GC\Marom\CranialNerves\Vagus\Analysis\';
fpath = 'C:\Users\devinomega\Dropbox\GC\Marom\CranialNerves\Vagus\Analysis\';
fname = 'demo_PT_pain.xlsx';

%%Path of Subject data
spath = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData';
sfiles = dir([spath '\TV*.mat']);

%%The Order of randomization
ordr = {'Active_25 Hz','Control_25 Hz', 'Active_100 Hz', 'Control_100 Hz', ...
    'Active_Burst', 'Control_Burst', 'Active_Bilateral 25 Hz', 'Control_Bilateral 25 Hz'};		

% Create a large file to save the data in Matlab
bigMat = [fpath 'tVNS_behavioral.mat'];
if exist(bigMat,'var') ~= 2
    % PT: (condition, subject,(PT/scalar/Dose))
    % NAS: (condition, subject, baseline rating) 
    tVNS_behavioral = struct('PT', zeros(8,30,3),'NAS',zeros(8,30,6));
    save(bigMat,'tVNS_behavioral');
else
    load(bigMat)
end

%create an anonymous function to change a letter of the alphabet to a number
Alphabet = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ';
let2num = @(c) 1+lower(c) - 'a';

% supress this warning that shows up later when we try and load a variable
% that's not there
warning('off', 'MATLAB:load:variableNotFound')
%%
% Create an OLE automation server for excel
% Specify file name
file = [fpath fname]; % This must be full path name
% Open Excel Automation server
Excel = actxserver('Excel.Application');
Workbooks = Excel.Workbooks;
% Make Excel visible
Excel.Visible=1;
% Open Excel file
Workbook=Workbooks.Open(file);
% Make the first sheet active
sheetnum=1;
Sheets = Excel.ActiveWorkBook.Sheets;
sheet1 = get(Sheets, 'Item', sheetnum);
invoke(sheet1, 'Activate');
Activesheet = Excel.Activesheet;

sheet2 = get(Sheets, 'Item', 2);
invoke(sheet2, 'Activate');
Activesheet2 = Excel.Activesheet;

% Grab the excel files where the data is stored
%   demo are the Demographics and the PT
%   NAS are the pain ratings
%   RND is the randomizations
[~, demoTxt, demoRaw] = xlsread([fpath fname],1);	
[~, NAStxt, NASraw] = xlsread([fpath fname],2);
[~, RNDtxt, RNDraw] = xlsread([fpath fname],3);	

[demoSubjR, demoSubjC] = find(strcmp(demoTxt,'Subject'));
[RNDSubjR, RNDSubjC] = find(strcmp(RNDtxt,'Subject'));
[NASSubjR, NASSubjC] = find(strcmp(NAStxt,'Subject'));
[~,demoCondC] =find(strcmp(demoTxt,ordr{1}));
[~,RNDCondC] =find(strcmp(RNDtxt,'Visit 1'));
[~,NASCondC] =find(strcmp(NAStxt,'Active_25 Hz'));
[~,condC] = cellfun(@(x) find(strcmp(demoTxt, x)),ordr,'UniformOutput', false);
%% Now go through all the subjects and create giant mat file and update
% excel
for i = 1:numel(sfiles)
    %extract the session and subject ID from the name
    curSubj = sfiles(i).name;
    session = regexp(curSubj,'TV_(\d\d)_(\d\d).mat','tokens','once');
    subjID = str2double(session{1});
    sessNum = str2double(session{2});   %The session 
    
    %Find the condition for that session # in the demo sheet and NAS
    conRNDCol = (sessNum-1) + RNDCondC; %find cond name in randomization sheet
    condLetter = RNDtxt{subjID+RNDSubjR, conRNDCol};  %condition as letter
    runNum = let2num(condLetter);
    cond = ordr{runNum};   % condition as string
    conDemoCol = (runNum-1)*3 + demoCondC;  %Find the run column in the demographics
    conNASCol = (runNum-1)*6 + NASCondC;    %Find the run column in the NAS
    conDemoRow = subjID+demoSubjR; %Find the run row for demo
    conNASRow = 1+subjID+NASSubjR; %Find the run row for demo
    
    % First - check to see if that subject already has values in their PT
    %   and NAS for that run
    if isnan(demoRaw{conDemoRow, conDemoCol})

        subjRun = load([spath '\' curSubj],'run');
        
        %Check to see if the condition in the .mat matches the
        %randomization
        if ~strcmp(subjRun.run,cond)
            uiwait(msgbox([curSubj '''s run and randomization do not match']))
            break
        end

        % load run_date, PT, final intensity, NAS 
        %   try to load resp - if it's there than the PT wasn't changed
        %   manually
        data = load([spath '\' curSubj],'run_Date','PT','pain','resp','run'); %load data

        %Handle the situaiton where the PT is a char
        if ischar(data.PT); data.PT = str2double(data.PT);end
        runInt = data.PT*2;
        
        %Check to see if the PT had changed (there will be no resp)
        if ~isfield(data, 'resp')   %check to see if the PT was changed
            data.PT = inputdlg(['PT changed. Enter PT manually for subj '  curSubj]);
        end
        
        if isempty(data.PT); break;end %case where dlg box is closed
        
        if iscell(data.PT); data.PT = str2double(data.PT);end %convert this to a double
        
        data.testNAS = inputdlg(['Enter test NAS for subject ' curSubj]); %Initial pain raiting
        
        if isempty(data.testNAS); break;end %case where dlg box is closed
        
        %calculate the values for the big file for each subject
        stimInt= [data.PT (runInt/data.PT) runInt]; %the PT and multiplier
        pain = [str2double(data.testNAS{:}) data.pain]; %all pain ratings
        sessNum = data.run; %Actual condition 
        runDate = data.run_Date; %Date the task was run

        % For each subject save their data into individual subject files
        % and put everything into a big file
        save([spath '\BIGMAT_' curSubj], 'stimInt','pain','sessNum')%,'-append')
        tVNS_behavioral.PT(runNum,subjID,:) = stimInt;
        tVNS_behavioral.NAS(runNum,subjID,:) = pain;
        
        % Specify sheet number, data, and range to write to for demo
        rangeDemo = [num2char(conDemoCol) num2str(conDemoRow) ':' num2char(conDemoCol+length(stimInt)-1) num2str(conDemoRow) ];
        rangeNAS = [num2char(conNASCol) num2str(conNASRow) ':' num2char(conNASCol+length(pain)-1) num2str(conNASRow) ];
        
        ActivesheetRange = get(Activesheet,'Range',rangeDemo);
        set(ActivesheetRange, 'Value', stimInt);
        ActivesheetRange = get(Activesheet2,'Range',rangeNAS);
        set(ActivesheetRange, 'Value', pain);
    end
end

% Turn the warning back on
warning('off', 'MATLAB:load:variableNotFound')
%% Save everything
save(bigMat, 'tVNS_behavioral')
% Save file
invoke(Workbook,'Save')
% Close Excel and clean up
invoke(Excel,'Quit');
delete(Excel);
clear Excel;