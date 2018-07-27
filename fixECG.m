% Certain subjects had the incorrect calculations for the BPM - this goes
% through those subjects and redoes the calculations
% Messed up HR average!

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
        if ~isfield(curData.ECG,'AvgBPM_5')
            curData.ECG.AvgBPM_5 = zeros(5,size(wndwsAvg_5,1)); % Avg BPM for [pre stim post] across a fixed window
        elseif size(curData.ECG.AvgBPM_5,2) ~= size(wndwsAvg_5,1)
            curData.ECG.AvgBPM_5 = zeros(5,size(wndwsAvg_5,1)); % Avg BPM for [pre stim post] across a fixed window
        end
        
        if ~isfield(curData.ECG,'medBPM_5')
            curData.ECG.medBPM_5 = zeros(5,size(wndwsAvg_5,1)); % Avg BPM for [pre stim post] across a fixed window
        end
        
        curData.ECG.AvgBPM_5 = zeros(5,size(wndwsAvg_5,1)); % Avg BPM for [pre stim post] across a fixed window
        
        for n = 1:numel(curData.ECG.locsData)
            locs = curData.ECG.locsData{n} ;
            if ~isnan(curData.ECG.locsData{n} )
                % ==================================================================
                % Average BPM across all periods (pre, stim,post)
                % ==================================================================

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
                
                curData.ECG.AvgBPM(n,:) = temp;    
                % ==================================================================
                % Average BPM across 5 sec bins
                % ==================================================================
%                 wndwsAvg = [[0:5:120]', [5:5:125]'];   %Time bins
%                 wndwSz = diff(wndwsAvg,1,2);    %length of time bins

                allLocs_5 = locs >= wndwsAvg_5(:,1) & locs < wndwsAvg_5(:,2); %Find the locs (peaks) within each bin
                splLocs = allLocs_5; %Added beat before and after
                edgeDetect = diff(allLocs_5,1,2); %detect edge
                splLocs(edgeDetect == 1) = 1;  %add beat immediately before
                splLocs([false(size(splLocs,1),1), edgeDetect] == -1) = 1; % Add the beat immediatly after

                temp = zeros(size(wndwsAvg_5,1),1);
                for t = 1:size(allLocs_5,1)
                    curLocs = locs(splLocs(t,:));
                    bop = (numel(curLocs)-3) + ((curLocs(2) - wndwsAvg_5(t,1))/(curLocs(2)-curLocs(1))) ...
                        + ((wndwsAvg_5(t,2)- curLocs(end-1))/(curLocs(end)-curLocs(end-1)));        %Average BPM - partial
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
                % ==================================================================
                % Max BPM for Stim
                % ==================================================================
                % Find the max BPM for the Stim time using a sliding window with 1
                % second overalp
                % time windows for stim and post (stops 10 seconds before end of
                % window
%                 ppWndws = {[15 65] [75 115]};
                for w = 1:numel(ppWndws)
%                     maxDiff = 0;
%                     minDiff = 300;
                    curWndw = ppWndws{w};
                    for t = curWndw(1):curWndw(2)
                        %Design the window (10s) and find the locations
                        window = [t t+10];
                        
                        Indx = find(locs >= window(1) &  locs <= window(2));
                        partLoc = locs([Indx(1)-1 Indx Indx(end)+1]);  % Add the beat immediatly before and after
                        part = (numel(Indx)-1) + (partLoc(2)-window(1))/(partLoc(2)-partLoc(1)) + ...
                            (window(2)-partLoc(end-1))/(partLoc(end)-partLoc(end-1));  %calculate partial beats in a period of time
                        
                        curBPM = (part/10)*60;
                        %                 locs(locs >= window(1) &  locs <= window(2));
                        %                 curLoc = mean(diff(locs(locs >= window(1) &  locs <= window(2))));
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
                    minBPM = minDiff - curData.ECG.AvgBPM(n,1) ;  %baseline correct
                    maxBPM = maxDiff - curData.ECG.AvgBPM(n,1);  %baseline correct
                    curData.ECG.minBPM(n,ind) = [minTime minBPM];
                    curData.ECG.maxBPM(n,ind) = [maxTime maxBPM];
                end
                
                allLocs_15 = locs >= wndwsAvg_15(:,1) & locs < wndwsAvg_15(:,2); %Find the locs (peaks) within each bin
                % ==================================================================
                % Heart Rate Variability
                % SDRR Standard deviation of RR intervals  
                % STD([RR1, RR2,...])
                % ==================================================================
                curData.ECG.hrVarData(n,:) = [std(diff(locs(allLocs(1,:)))) std(diff(locs(allLocs(2,:))))...
                    std(diff(locs(allLocs(3,:))))];
                
                % Fifteen seconds windows
                temp = zeros(size(wndwsAvg_15,1),1);  
                for t = 1:size(allLocs_15,1)
                    temp(t) = std(diff(locs(allLocs_15(t,:))));
                end
                curData.ECG.hrVarData_15(n,:) = temp;  
                % ==================================================================
                % RMSSD - square root of the mean squared difference of
                % successive NN intervals
                % sqrt(mean([(RR1-RR2)^2,(RR2-RR3)^2...]))
                % ==================================================================
            
                curData.ECG.RMSSD(n,:) = [sqrt(mean(diff(diff(locs(allLocs(1,:)))).^2)) ...
                    sqrt(mean(diff(diff(locs(allLocs(2,:)))).^2)) ...
                    sqrt(mean(diff(diff(locs(allLocs(3,:)))).^2))];
                
                % Fifteen seconds windows
                temp = zeros(size(wndwsAvg_15,1),1);  
                for t = 1:size(allLocs_15,1)
                    temp(t) = sqrt(mean(diff(diff(locs(allLocs_15(t,:)))).^2));
                end
                curData.ECG.RMSSD_15(n,:) = temp;  
                % ==================================================================
                % pRR50 percent of succesive RR intervals greater than 50ms
                % ==================================================================
                curData.ECG.pRR50(n,:) = [sum(diff(diff(locs(allLocs(1,:))))>.05)/numel(diff(diff(locs(allLocs(1,:))))) ...
                    sum(diff(diff(locs(allLocs(2,:))))>.05)/numel(diff(diff(locs(allLocs(2,:))))) ...
                    sum(diff(diff(locs(allLocs(3,:))))>.05)/numel(diff(diff(locs(allLocs(3,:)))))];
                
                % Fifteen seconds windows
                temp = zeros(size(wndwsAvg_15,1),1);  
                for t = 1:size(allLocs_15,1)
                    temp(t) = sum(diff(diff(locs(allLocs_15(t,:))))>.05)/numel(diff(diff(locs(allLocs(1,:)))));
                end
                curData.ECG.pRR50_15(n,:) = temp;  
            else
                curData.ECG.AvgBPM(n,:) = locs;
                curData.ECG.AvgBPM_5(n,:) = locs;
                curData.ECG.medBPM(n,:) = locs;
                curData.ECG.minBPM(n,:) = locs;
                curData.ECG.maxBPM(n,:) = locs;
                curData.ECG.hrVarData(n,:) = locs;
                curData.ECG.RMSSD(n,:)= locs;
                curData.ECG.pRR50(n,:)= locs;
                curData.ECG.hrVarData_15(n,:) = locs;
                curData.ECG.RMSSD_15(n,:)= locs;
                curData.ECG.pRR50_15(n,:)= locs;
            end
        end
        %     Only save the ECG data if all blocks are done - otherwise you'll have
        %     partially complete data
        ECG = curData.ECG;
        save([fPath curSubj], 'ECG', '-append')
        clear ECG
    end
end
