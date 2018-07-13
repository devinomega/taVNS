%% ECG Processing Pipeline
%
% 1) Update the tables to include pain scale ratings and thresholds
% 2) Change the ECG files names (.cnt & .evt) to be in run format 
% 3) Peak Pick

%% PreProcessing
updateTables

updateECG

%% ECG
%heart rate picking for the five blocks before during and after stim
ECG_preProc('skip',false,'startName','TV_20_02')
% heart rate picking for the baseline 2.5 min
ECGbaselinePreProc('skip',false,'startName','TV_27_04')
%heart rate picking for the five blocks post number rating
ECGpostPreProc('skip',false,'startName','TV_02_01')
%% Resp
% resperation picking for the five blocks before during and after stim
breath_preProc('skip',false,'startName','TV_20_02')
% resp rate picking for the baseline 2.5 min
breathBaselinePreProc('skip',false,'startName','TV_27_08')
%% Run Stats

HR_Stats()