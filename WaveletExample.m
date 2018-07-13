%% Load and plot an ECG waveform where the R peaks of the QRS complex have been annotated by two or more cardiologists. The ECG data and annotations are taken from the MIT-BIH Arrythmia Database. The data are sampled at 360 Hz.
%Segmenting
dataLngth = 135*2000;   %time*sampling rate 

% parameters
sampleRate = 2000; % Hz
lowEnd = .5; % Hz
highEnd = 26; % Hz
tm = 0:(1/2000):135;
tm(1) = [];

%load data
data = pop_loadeep_v4('C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData\TV_03_01.cnt');   %load their eeg data
curData = double(data.data(33,:)); %Raw ECG
pre = [data.event(~cellfun(@isempty,regexp({data.event.type}, ...
            '1\d{1}1'))).latency];
ECG.allData = cellfun(@(i) curData(i:(i+dataLngth-1)), num2cell(pre),'UniformOutput',false);

temp = ECG.allData{1};
ecgsig_real= temp(1:270000);   
% ecgSig_inv = ecgsig_real*-1;    

filterOrder = 2; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
[b, a] = butter(filterOrder, [lowEnd highEnd]/(sampleRate/2)); % Generate filter coefficients
filtData = filtfilt(b, a,ecgsig_real); % Apply filter to data using zero-phase filtering

% filterOrder = 2; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
% [b, a] = butter(filterOrder, [lowEnd highEnd]/(sampleRate/2)); % Generate filter coefficients
% filtData = filtfilt(b, a,ecgSig_inv); % Apply filter to data using zero-phase filtering

figure
fig1  = plot(tm(1:10000),filtData(1:10000));   % 20 seconds of data
% hold on
% plot(tm(ann),ecgsig(ann),'ro')
xlabel('Seconds')
ylabel('Amplitude')
title('Subject - TV-03-01')

%% You can use wavelets to build an automatic QRS detector for use in applications like R-R interval estimation.
% There are two keys for using wavelets as general feature detectors:
% The wavelet transform separates signal components into different frequency bands enabling a sparser representation of the signal.
% You can often find a wavelet which resembles the feature you are trying to detect.
% The 'sym4' wavelet resembles the QRS complex, which makes it a good choice for QRS detection. To illustrate this more clearly, extract a QRS complex and plot the result with a dilated and translated 'sym4' wavelet for comparison.
qrsEx = ecgsig_real(17800:19600);
qrsEx = qrsEx-mean(qrsEx(1:680)); %baseline correct  move it on up 0
% norm_data = (max(blc)-ymin(blc))*(x(-xmin)/(xmax-xmin) + ymin);
[mpdict,~,~,longs] = wmpdictionary(numel(qrsEx),'lstcpt',{{'sym4',5}}); %level 3, wavelet sym4
figure
plot(qrsEx)
hold on
plot(500*circshift(mpdict(:,11),[800 0]),'r') %11 is the number of 
axis tight
legend('QRS Complex','Sym4 Wavelet')
title('Comparison of Sym4 Wavelet and QRS Complex')

%% Use the maximal overlap discrete wavelet transform (MODWT) to enhance the R peaks in the ECG waveform. The MODWT is an undecimated wavelet transform, which handles arbitrary sample sizes.
% First, decompose the ECG waveform down to level 5 using the default 'sym4' wavelet. Then, reconstruct a frequency-localized version of the ECG waveform using only the wavelet coefficients at scales 4 and 5. The scales correspond to the following approximate frequency bands.

% Scale 4 -- [11.25, 22.5) Hz
% Scale 5 -- [5.625, 11.25) Hz.

%  This covers the passband shown to maximize QRS energy.
wt = modwt(filtData,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');

figure
for r = 1:(size(wt,1)-1)
   subplot(5,1,r) 
   plot(wt(1,:))
end

wt = modwt(filtData,'sym6');
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');
%% Use the squared absolute values of the signal approximation built from the wavelet coefficients and employ a peak finding algorithm to identify the R peaks.

y = abs(y).^2;
[qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',50,...
    'MinPeakDistance',0.150);
figure
plot(tm,y)
hold on
plot(locs,qrspeaks,'ro')
xlabel('Seconds')
title('R Peaks Localized by Wavelet Transform with Automatic Annotations')


figure
plot(tm,filtData);   % 20 seconds of data
hold on
plot(locs,qrspeaks,'ro')
xlabel('Seconds')
title('R Peaks Localized by Wavelet Transform with Automatic Annotations')
%% At the command line, you can compare the values of tm(ann) and locs, which are the expert times and automatic peak detection times respectively. Enhancing the R peaks with the wavelet transform results in a hit rate of 100% and no false positives. The calculated heart rate using the wavelet transform is 88.60 beats/minute compared to 88.72 beats/minute for the annotated waveform.
% If you try to work on the square magnitudes of the original data, you find the capability of the wavelet transform to isolate the R peaks makes the detection problem much easier. Working on the raw data can cause misidentifications such as when the squared S-wave peak exceeds the R-wave peak around 10.4 seconds.

figure
plot(tm,ecgsig_real,'k--')
hold on
plot(tm,y,'r','linewidth',1.5)
plot(tm,abs(ecgsig_real).^2,'b')
plot(tm(ann),ecgsig_real(ann),'ro','markerfacecolor',[1 0 0])
set(gca,'xlim',[10.2 12])
legend('Raw Data','Wavelet Reconstruction','Raw Data Squared', ...
    'Location','SouthEast');
xlabel('Seconds')

%% Using findpeaks on the squared magnitudes of the raw data results in twelve false positives.
[qrspeaks,locs] = findpeaks(ecgsig_real.^2,tm,'MinPeakHeight',0.35,...
    'MinPeakDistance',0.150);

%% In addition to switches in polarity of the R peaks, the ECG is often corrupted by noise.
load mit200
figure
plot(tm,ecgsig)
hold on
plot(tm(ann),ecgsig(ann),'ro')
xlabel('Seconds')
ylabel('Amplitude')
title('Subject - MIT-BIH 203 with Expert Annotations')

wt = modwt(ecgsig,5);
wtrec = zeros(size(wt));
wtrec(4:5,:) = wt(4:5,:);
y = imodwt(wtrec,'sym4');
y = abs(y).^2;
[qrspeaks,locs] = findpeaks(y,tm,'MinPeakHeight',0.1,...
    'MinPeakDistance',0.150);

%% Use the MODWT to isolate the R peaks. Use findpeaks to determine the peak locations. Plot the R-peak waveform along with the expert and automatic annotations.

figure
plot(tm,y)
title('R-Waves Localized by Wavelet Transform')
hold on
hwav = plot(locs,qrspeaks,'ro');
hexp = plot(tm(ann),y(ann),'k*');
xlabel('Seconds')
legend([hwav hexp],'Automatic','Expert','Location','NorthEast');

