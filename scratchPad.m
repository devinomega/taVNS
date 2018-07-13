%% Testing some stuff out
% Scratch Pad

curSubj = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData\TV_10_01.cnt'; %subject currently being worked on    
data = pop_loadeep_v4(curSubj);   %load their eeg data
   
%For pre, stim, post, and end of post determine the index where each
%block's trigger is located along all the data
pre = [data.event(~cellfun(@isempty,regexp({data.event.type}, ...
    '1\d{1}1'))).latency];
stim = [data.event(~cellfun(@isempty,regexp({data.event.type}, ...
    '1\d{1}2'))).latency];
post = [data.event(~cellfun(@isempty,regexp({data.event.type}, ...
    '1\d{1}3'))).latency];
postEnd = [data.event(~cellfun(@isempty,regexp({data.event.type}, ...
    '1\d{1}4'))).latency];

curData = double(data.data(33,:)); %Raw ECG

ECG.allData = cellfun(@(i,j) curData(i:j), num2cell(pre), ...
    num2cell(postEnd),'UniformOutput',false ); %Data not interpolate

tmLn = 2000*135;

tmp = ECG.allData{1};
tmp = tmp(1:tmLn);
rawTime  = (0:numel(tmp)-1)/2000;

figure
plot(tmp)


wt = modwt(tmp,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');
MODWT = abs(y).^2;
threshold = quantile(MODWT,.98);
[peaks,locs] = findpeaks(MODWT, rawTime,'MinPeakHeight',...
    threshold,'MinPeakDistance',0.5);

figure
plot(rawTime, MODWT)
title('R-Waves Localized by Wavelet Transform')
hold on
plot(locs,peaks,'ro')
xlabel('Seconds')

y = fft(tmp);
y(1) = [];
figure
plot(y,'ro')
xlabel('real(y)')
ylabel('imag(y)')
title('Fourier Coefficients')
n = length(y);
power = abs(y(1:floor(n/2))).^2; % power of first half of transform data
% maxfreq = 1/2;                   % maximum frequency
% freq = (1:n/2)/(n/2); %*maxfreq;    % equally spaced frequency grid
 freq = 2000*(1:(n/2))/n;
 figure
plot(freq,power)
xlabel('Hz')
ylabel('Power')

sampleRate = 2000; % Hz
lowEnd = .5; % Hz
highEnd = 30; % Hz
filterOrder = 2; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
[b, a] = butter(filterOrder, [lowEnd highEnd]/(sampleRate/2)); % Generate filter coefficients
filteredData = filtfilt(b, a, tmp); % Apply filter to data using zero-phase filtering
figure
plot(filteredData)



wt = modwt(filteredData,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');
MODWT = abs(y).^2;
threshold = quantile(MODWT,.98);
[peaks,locs] = findpeaks(MODWT, rawTime,'MinPeakHeight',...
    threshold,'MinPeakDistance',0.5);

figure
plot(rawTime, MODWT)
title('R-Waves Localized by Wavelet Transform')
hold on
plot(locs,peaks,'ro')
xlabel('Seconds')


handles.threshold = quantile(handles.rawData,.975);
[handles.peaks,handles.locs] = findpeaks(handles.rawData,...
    handles.rawTime,'MinPeakHeight',handles.threshold,'MinPeakDistance',0.5);
    
y = fft(ECG.interpData(n,:));
y(1) = [];
plot(y,'ro')
xlabel('real(y)')
ylabel('imag(y)')
title('Fourier Coefficients')

n = length(y);
power = abs(y(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
plot(freq,power)
xlabel('Cycles/Year')
ylabel('Power')

%Find the Bpm of data
locs =peakPick(ECG.interpData(n,:),[curSubj(1:end-4) '_B' num2str(n)]);

% curLen =  numel(ECGallData{5});
% xx = (0:(135*2000))/2000;
% x = (0:(curLen-1))/2000;
% a = linspace(0,curLen,135*2000)/2000;
% test = zeros(1,curLen);
% test(stim(5) - pre(5)+1) = 1;
% test(post(5) - pre(5)+1) = 1;
% y = test;
% new = interp1(x,y,a);

%     set(f,'visible','on')

    % 
%     if numel(data.event)  < 26
%         answer = questdlg(['The triggers are weird for ' curSubj],...
%             'Trigger error','OK, continue','Break','OK, continue');
%         if strcmp(answer, 'Break')
%             return
%         end
%     end

%         preLoc = {RR_Peak_Pick(curBlock(5*2000+1:2000*15))};
%         stimLoc = {RR_Peak_Pick(curBlock(15*2000+1:2000*75))};
%         postLoc =  {RR_Peak_Pick(curBlock(75*2000+1:end))};
%         bpmData(n) = {60./diff(locs)}; %1./diff(locs(2:end-1)); 
        
        %graphing/fitting curve
%         xData = locs(2:end);
%         yData = bpm;
%         yy = spline(xData,yData,xx);

%         stimTrig = [xx(stim(n)-pre(n)+1) xx(post(n)-pre(n)+1)];
%         baseline = [(stim(n)-pre(n)+1-2000*10),(stim(n)-pre(n)+1)];
% %         baseCorrect = yy-mean(yy(baseline));
%         baseData(n,:)  = yy-mean(yy(baseline));
        
        %baseCorrect(baseline(1):end);
        
%         subplot(2,3,n);
%         plot(xx,yy)
%         xlabel('Seconds')
%         ylabel('BPM')
%         title(['Block ' num2str(n)])
%         vline(stimTrig,{'r','r'},{'Stim Start','Stim End'})
% %
% % DC Offset
% dc_off = ecgSig - mean(ecgSig);
% 
% % Detrend
% % find the mean
% detrend_data = detrend(ecgSig);
% figure
% plot(tmSec,ecgSig);
% hold on
% plot(tmSec,detrend_data)
% 
% % Filter pass
% 
% % Plot the data
% figure
% hax = axes;
% plot(tmSec(1:20000),ecgSig(1:20000))
% hold on
% yLimAx = get(hax,'YLim');
% 
% for i = 1:numel(trg)
%     plot(tmSec(trg(i)),yLimAx)
% end
% 
% %         curBlock = curBlock-mean(curBlock);
% %         curBlock = detrend(curBlock);
% % ECG analysis
% ecgSigEx = double(ecgSig(1:20000));
% tm2 = tmSec(1:20000);
% 
% 
% % Use the MODWT to isolate the R peaks 
% %   Use findpeaks to determine the peak locations. 
% %   Plot the R-peak waveform along with the expert and automatic annotations
% wt = modwt(ecgSigEx,5);
% wtrec = zeros(size(wt));
% wtrec(4:5,:) = wt(4:5,:);
% y = imodwt(wtrec,'sym4');
% y = abs(y).^2;
% [qrspeaks,locs] = findpeaks(y,tm2,'MinPeakHeight',.001,'MinPeakDistance',0.5);
% figure
% plot(tm2,y)
% title('R-Waves Localized by Wavelet Transform')
% hold on
% hwav = plot(locs,qrspeaks,'ro');
% xlabel('Seconds')
% 
% % %Detrend
% % detrend_data = detrend(ecgSigEx);
% % 
% % % DC Offset
% % dc_off = detrend_data- mean(detrend_data);
% % 
% % %normalize
% % norm_data = (dc_off - min(dc_off)) / ( max(dc_off) - min(dc_off) );
% % 
% % %filter
% % sampleRate = 2000; % Hz
% % lowEnd = .1; % Hz
% % highEnd = 40; % Hz
% % filterOrder = 2; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
% % [b, a] = butter(filterOrder, [lowEnd highEnd]/(sampleRate/2)); % Generate filter coefficients
% % filteredData = filtfilt(b, a, norm_data); % Apply filter to data using zero-phase filtering
% % FreqzA=[0.2 4 8 15]/1000;
% % FreqzB=[4 8 12 30]/1000;
% % [b,a]=butter(2,[FreqzA(gg) FreqzB(gg)],'bandpass');
% % or 
% % [b,a]=butter(2,[0.5 8]/1000,'bandpass')
% % figure
% % plot(tm2,norm_data)
% % hold on
% % plot(tm2,filteredData)
% 
% %Another option where you use the shape of the QRS wave to fit it to the
% %data
% 
% % [mpdict,~,~,longs] = wmpdictionary(numel(ecgSigEx),'lstcpt',{{'sym4',3}});
% % figure
% % plot(ecgSigEx)
% % hold on
% % plot(2*circshift(mpdict(:,11),[-2 0]),'r')
% % axis tight
% % legend('QRS Complex','Sym4 Wavelet')
% % title('Comparison of Sym4 Wavelet and QRS Complex')
% % 
% % plot(tmSec(trg(1:26)), ecgSig(trg(1:26)),'ro')



% Make a table of actual data
% Make a table of the measures
% ==================================================================
% GroupXBlockXPrePost
% ==================================================================

% curMeas = cell(numel(ordr)*numel(blck)*numel(prePost),1);  %variable names
% curWiMat = zeros(numel(curMeas),3);
% 
% m=1;
% for o = 1:numel(ordr)
%     for b = 1:numel(blck)
%         for p = 1:numel(prePost)
%             curWiMat(m,:) = [o, b, p];
%             curMeas{m} = [ordr{o} '_0' blck{b} '_' prePost{p}];
%             m = m+1;
%         end
%     end
% end
% curMeas = strrep(curMeas,' ','_'); 
% % dataShape = permute(allData,[2 4 3 1]);
% dataShape = reshape(dataShape,[numSubj numel(curMeas)]);
% dataShape = dataShape(~all(isnan(dataShape),2),:);
% dataShape(isnan(dataShape)) = 999;
% 
% xlswrite([sPath 'Group_Block_PreStim.xlsx'],[curMeas';num2cell(dataShape)]);
% % ==================================================================
% % Difference between Pre and Stim
% % ==================================================================
% curMeas = cell(numel(ordr)*numel(blck),1);  %variable names
% curWiMat = zeros(numel(curMeas),2);
% 
% m=1;
% for o = 1:numel(ordr)
%     for b = 1:numel(blck)
%             curWiMat(m,:) = [o, b];
%             curMeas{m} = [ordr{o} '_0' blck{b} '_diff'];
%             m = m+1;
%     end
% end
% 
% diffData = diff(allData,1,4);
% dataShape = permute(diffData,[2 4 3 1]);
% dataShape = reshape(dataShape,[numSubj numel(curMeas)]);
% dataShape = dataShape(~all(isnan(dataShape),2),:);
% dataShape(isnan(dataShape)) = 999;
% 
% xlswrite([sPath  'Group_Block_Diff.xlsx'],[curMeas';num2cell(dataShape)]);
% 
% % 9999
% % xlswrite([sPath  'HRV_Group_Block_Diff.xlsx'],[meas';num2cell(dataShape)]);
% xlswrite([sPath  'RMSSD_Group_Block_Diff.xlsx'],[curMeas';num2cell(dataShape)]);
% % xlswrite([sPath  'pRR50_Group_Block_Diff.xlsx'],[meas';num2cell(dataShape)]);
% 
% % ==================================================================
% % Difference between Pre and Stim collapse block
% % ==================================================================
% curMeas = ordr;  %variable names
% curWiMat = 1:numel(ordr);
% 
% diffData = diff(allData,1,4);
% % dataShape = permute(nanmean(diffData,3),[2 1]);
% dataShape = dataShape(~all(isnan(dataShape),2),:);
% dataShape(isnan(dataShape)) = 999;
% 
% xlswrite([sPath  'Group_Diff.xlsx'],[curMeas';num2cell(dataShape)]);
% 
% xlswrite([sPath  'HRV_Group_Diff.xlsx'],[curMeas;num2cell(dataShape)]);
% 
% % ==================================================================
% % 5 Second bins
% % ==================================================================
% curMeas = cell(numel(ordr)*numel(blck)*(numel(prePost)-1),1);  %variable names
% curWiMat = zeros(numel(curMeas),3);
% 
% m=1;
% for o = 1:numel(ordr)
%     for b = 1:numel(blck)
%         for p = 2:numel(prePost)
%             curWiMat(m,:) = [o, b, p];
%             curMeas{m} = [ordr{o} '_0' blck{b} '_' prePost{p}];
%             m = m+1;
%         end
%     end
% end
% curMeas = strrep(curMeas,' ','_'); 
% baseChunk = (allData(:,:,:,2:13))-repmat(allData(:,:,:,1),1,1,1,12);
% % dataShape = permute(baseChunk,[2 4 3 1]);
% dataShape = reshape(dataShape,[numSubj numel(curMeas)]);
% dataShape = dataShape(~all(isnan(dataShape),2),:);
% dataShape(isnan(dataShape)) = 999;
% 
% xlswrite([sPath 'Group_Block_5S_Bins.xlsx'],[curMeas';num2cell(dataShape)]);
% 
% % ==================================================================
% % 5 Second bins no block
% % ==================================================================
% curMeas = cell(numel(ordr)*(numel(prePost)-1),1);  %variable names
% curWiMat = zeros(numel(curMeas),2);
% 
% m=1;
% for o = 1:numel(ordr)
%         for p = 2:numel(prePost)
%             curWiMat(m,:) = [o, p];
%             curMeas{m} = [ordr{o} '_' prePost{p}];
%             m = m+1;
%         end
% end
% curMeas = strrep(curMeas,' ','_'); 
% baseChunk = (allData(:,:,:,2:13))-repmat(allData(:,:,:,1),1,1,1,12);
% subjGrpTime = squeeze(nanmean(baseChunk,3));
% % dataShape = permute(subjGrpTime,[2 3 1]);
% dataShape = reshape(dataShape,[numSubj numel(curMeas)]);
% dataShape = dataShape(~all(isnan(dataShape),2),:);
% dataShape(isnan(dataShape)) = 999;
% 
% xlswrite([sPath 'Group_5S_Bins.xlsx'],[curMeas';num2cell(dataShape)]);
%% Bin
% Bin the heart rates into 5/10 seconds bins
% 
% wndwsAvg = [[5, 15:5:120]', [15, 20:5:125]'];   %Time bins
% wndwSz = diff(wndwsAvg,1,2);    %length of time bins
% z = 1;
% for i = z:numel(fFiles)
%     curSubj = fFiles(i).name;
%     load([fPath curSubj],'sessNum','ECG')   %Load up the ECG file
%     ECG.AvgBPM_5 = zeros(5,size(wndwsAvg,1)); %Add in the 5 second bin
%     
%     %Go through each block
%     for b = 1:numel(ECG.locsData)
%         locs = ECG.locsData{b}; %current location data
%         if ~isnan(ECG.locsData{b})
%             allLocs = locs >= wndwsAvg(:,1) & locs < wndwsAvg(:,2); %Find the locs (peaks) within each bin
%             splLocs = allLocs; %Added beat before and after
%             edgeDetect = diff(allLocs,1,2); %detect edge
%             splLocs(edgeDetect == 1) = 1;  %add beat immediately before
%             splLocs([false(size(splLocs,1),1), edgeDetect] == -1) = 1; % Add the beat immediatly after
%             
%             temp = zeros(size(wndwsAvg,1),1);
%             for n = 1:size(allLocs,1)
%                 curData = locs(splLocs(n,:));
%                 bop = (numel(curData)-3) + ((curData(2) - wndwsAvg(n,1))/(curData(2)-curData(1))) ...
%                     + ((wndwsAvg(n,2)- curData(end-1))/(curData(end)-curData(end-1)));        %Average BPM - partial
%                 temp(n) = (bop/wndwSz(n))*60;
%             end
%             ECG.AvgBPM_5(b,:) = temp;    
%         else
%             ECG.AvgBPM_5(b,:) = locs;
%         end
%     end
%     save([fPath curSubj], 'ECG', '-append')
% end
