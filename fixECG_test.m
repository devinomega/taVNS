% There are several ways to calculate BPM in a fixed windo (e.g. 5 sec)
% This is based on Chabot (1991) & Graham (1978)
%
%   1) Count number of beats (including fractions of beats w/ a time period
%   2) Counting QRS complexes (number of R's no partial beats
%   3) Average all IBIs in a period including weighted fractions (???)
%   4) Calculate based on time period for beats w/in fixed time window
%   5) 

%
fPath = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\taVNS_SubjData071118\';
% fPath='C:\Users\Dadair\Dropbox\MATLAB\tVNS\NI\subjData\';
fFiles = dir([fPath '\BIGMAT*.mat']);
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



for i = 1:numel(fFiles)
    curSubj = fFiles(i).name; %subject currently being worked on
    curData = load([fPath curSubj],'ECG');
    
    %     Check to see if the subject has an ECG variable
    if isfield(curData, 'ECG')
        locs = curData.ECG.locsData{n} ;
        for n = 1:numel(curData.ECG.locsData)
            % ==================================================================
            % Average BPM across 5 sec bins
            % ==================================================================
            %This does not work for the first or last time segment
            allLocs_5 = locs >= wndwsAvg_5(:,1) & locs < wndwsAvg_5(:,2); %Find the locs (peaks) within each bin
            splLocs = allLocs_5; %Added beat before and after
            edgeDetect = diff(allLocs_5,1,2); %detect edge
            splLocs(edgeDetect == 1) = 1;  %add beat immediately before
            splLocs([false(size(splLocs,1),1), edgeDetect] == -1) = 1; % Add the beat immediatly after
            
            temp = zeros(size(wndwsAvg_5,1),1);
            for t = 1:size(allLocs_5,1)
                bufLocs = locs(splLocs(t,:));   %locs with with the previous and next beat added
                timeConst_Locs = locs(allLocs_5(t,:));  %Beats w/in the time frame
                if timeConst_Locs <
                bop = (numel(bufLocs)-3) + ((bufLocs(2) - wndwsAvg_5(t,1))/(bufLocs(2)-bufLocs(1))) ...
                    + ((wndwsAvg_5(t,2)- bufLocs(end-1))/(bufLocs(end)-bufLocs(end-1)));        %Average BPM - partial
                temp(t) = (bop/wndwSz_5(t))*60;
            end
            
            curData.ECG.AvgBPM_5(n,:) = temp;
            % ==================================================================
            % Median BPM across 5 sec bins
            % ==================================================================
            %                 wndwsAvg = [[0:5:120]', [5:5:125]'];   %Time bins
            %                 wndwSz = diff(wndwsAvg,1,2);    %length of time bins
            
            allLocs_5 = locs >= wndwsAvg_5(:,1) & locs < wndwsAvg_5(:,2); %Find the locs (peaks) within each bin
            temp = zeros(size(wndwsAvg_5,1),1);
            
            for t = 1:size(allLocs_5,1)
                temp(t) = 60./median(diff(locs(allLocs_5(t,:))));
            end
            
            curData.ECG.medBPM_5(n,:) = temp;
        end
    end
end
