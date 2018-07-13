%% ECG Processing Pipeline
%
% 1) Update the tables to include pain scale ratings and thresholds
% 2) Change the ECG files names (.cnt & .evt) to be in run format 
% 3) Peak Pick

%% PreProcessing
%
updateTables

updateECG

ECG_preProc('filePath', 'C:\Users\bikson\Documents\taVNS_data\subjData\', 'skip',false,'startName','TV_11_07')
%% heart rate picking for the five blocks^
breath_preProc('filePath', 'C:\Users\bikson\Documents\taVNS_data\subjData\', 'skip',false)
%%resperation picking for the five blocks^
ECGbaselinePreProc('filePath', 'C:\Users\bikson\Documents\taVNS_data\subjData\', 'skip',false,'startName','TV_27_04')
%% heart rate for baseline^
breathBaselinePreProc('filePath', 'C:\Users\bikson\Documents\taVNS_data\subjData\', 'skip',false,'startName','TV_27_08')
%% Run Stats

Hr_Stats()
