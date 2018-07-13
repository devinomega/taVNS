% Redo the calculation of heart rate - was done incorrectly the first time

clear
curSubj = 'BIGMAT_TV_02_02.mat';
load(curSubj,'ECG')


ECG.AvgBPM = zeros(5,3); %Avg BPM for [pre stim post]
ECG.minBPM = zeros(5,6);  %min BPM in a sliding mindow for stim and post
ECG.minBPM = zeros(5,6);   % same with min (first 3 are stim, second 3 are post)

for n = 1:5
    locs = ECG.locsData{n};
    
    %Pre, stim and Post segments of peak locations
    preLoc = locs(locs >= 5 & locs < 15);
    stimLoc = locs(locs >= 15 & locs < 75);
    postLoc = locs(locs >= 75);
    
    ECG.AvgBPM(n,:) = [mean(diff(preLoc)*60) mean(diff(stimLoc)*60) ...
        mean(diff(postLoc)*60)];
    
    % time windows for stim and post
    ppWndws = {[15 65] [65 125]};
    for w = 1:numel(ppWndws)
        maxDiff = 0;
        minDiff = 300;
        curWndw = ppWndws{w};
        for t = curWndw(1):curWndw(2)
            %Design the window (10s) and find the locations
            window = [t t+10];
            locs(locs >= window(1) &  locs <= window(2));
            curLoc = mean(diff(locs(locs >= window(1) &  locs <= window(2))));
            %Update current max
            if curLoc > maxDiff
                maxTime = window;
                maxDiff = curLoc;
            end
            if curLoc < minDiff
                minTime = window;
                minDiff = curLoc;
            end
        end
        indx = (w-1)*3+1:(w-1)*3+3;
        minBPM = minDiff*60 - ECG.AvgBPM(n,1) ;  %baseline correct
        maxBPM = maxDiff*60 - ECG.AvgBPM(n,1);  %baseline correct
        ECG.minBPM(n,indx) = [minTime minBPM];
        ECG.maxBPM(n,indx) = [maxTime maxBPM];
    end
end

save(curSubj,'ECG','-append')