%% taVNS Eyelink Analysis Script
%
% CNT file - physiological recording from ANT
% EDF file - eyetracking data
% MAT file - basic information form the run

fPath = 'C:\Users\dadair\Dropbox\MATLAB\tVNS\NI\subjData\';
subjID = 'TV_01_01'; 

edf0 = Edf2Mat( [subjID '.edf']);
areaC = strcmp('AREA', edf0.Events.pupilInfo{:});

%
stimTime = edf0.Events.Messages.time(~cellfun(@isempty,regexp(edf0.Events.Messages.info,'STIM*')));
stimNorm = stimTime-(edf0.Samples.time(1)+1);

normT = edf0.normalizedTimeline/1000;
% xLim = linspace(normT(1),round(normT(end),-3),11)/1000;

%
fig = figure;
hax = axes;
x = (0:length(edf0.Samples.pa))/(500*60);
x(end)=[];
data = edf0.Samples.pa(:,1) ;

plot(x,edf0.Samples.pa(:,1))
% plot(edf0.Samples.pa(:,1))
hold on
yLimAx = get(hax,'YLim');
xLimAx= get(hax,'XLim');
colors = {'r','g'}; %{end, beg}
xlim([0 x(end)])
xlabel('Time (minutes)')
ylabel('Pupil area (arbitrary)')

dataN = stimNorm/(1000*60);
for i = 1:length(stimNorm)
%     plot([stimTime(i) stimTime(i)],yLimAx, colors{mod(i,2)+1},'LineWidth',1.5)
%      plot([stimNorm(i)*2 stimNorm(i)*2],yLimAx, colors{mod(i,2)+1},'LineWidth',1.5)
   plot([dataN(i) dataN(i)],yLimAx, colors{mod(i,2)+1},'LineWidth',1.5)
end

%%
base = edf0.Samples.pa(stimNorm(1)-999:stimNorm(1),1);
stim =  edf0.Samples.pa(stimNorm(1):stimNorm(2),1);

data(data==0) = NaN;

figure
plot(x,data,':')
% plot(edf0.Samples.pa(:,1))
hold on
yLimAx = get(hax,'YLim');
xLimAx= get(hax,'XLim');
colors = {'r','g'}; %{end, beg}
xlim([0 x(end)])
xlabel('Time (minutes)')
ylabel('Pupil area (arbitrary)')

dataN = stimNorm/(1000*60);
for i = 1:length(stimNorm)
%     plot([stimTime(i) stimTime(i)],yLimAx, colors{mod(i,2)+1},'LineWidth',1.5)
%      plot([stimNorm(i)*2 stimNorm(i)*2],yLimAx, colors{mod(i,2)+1},'LineWidth',1.5)
   plot([dataN(i) dataN(i)],yLimAx, colors{mod(i,2)+1},'LineWidth',1.5)
end

%% Pupilometry analysis

% setup path
% clear; clc; close all;
% thispath = '~/Dropbox/code/PupilPreprocessing/';
% addpath(thispath);

% first, define the asc filename
edfFile = 'TV_02_02.edf.edf';
ascFile = regexprep(edfFile, 'edf', 'asc');

% set path to FieldTrip - get this from http://www.fieldtriptoolbox.org/download
%addpath('~/Documents/fieldtrip/');
ft_defaults;

% ============================================== %
% 1. convert edf file to asc
% ============================================== %
% 
% Have to do this in BASH
%
% if ~exist(ascFile, 'file'),
%     if ismac,
%         edf2ascPath = [thispath '/edf2asc-mac'];
%     elseif isunix,
%         edf2ascPath = [thispath '/edf2asc-linux'];
%     else
%         error('Sorry, I don''t have an edf2asc converter for Windows')
%     end
%     % use this converter to create the asc file
%     % failsafe mode avoids errors when some samples are missing
%     system(sprintf('%s %s -failsafe -input', edf2ascPath, edfFile));
% end
% assert(exist(ascFile, 'file') > 1, 'Edf not properly converted');

% ============================================== %
% 2. create a FieldTrip-style data structure
% ============================================== %

% read in the asc EyeLink file
asc = read_eyelink_ascNK_AU(ascFile);

% create events and data structure, parse asc
[data, event, blinksmp, saccsmp] = asc2dat(asc);

% ============================================== %
% 3. interpolate Eyelink-defined and additionally detected blinks
% ============================================== %

plotMe = true;
newpupil = blink_interpolate(data, blinksmp, plotMe);
data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:) = newpupil;

% ============================================== %
% (optional): regress out blink- and saccade-linked pupil response
% http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0155574 
% ============================================== %

pupildata = data.trial{1}(~cellfun(@isempty, strfind(lower(data.label), 'eyepupil')),:);
newpupil = blink_regressout(pupildata, data.fsample, blinksmp, saccsmp, 1, 1);
% put back in fieldtrip format
data.trial{1}(~cellfun(@isempty, strfind(lower(data.label), 'eyepupil')),:) = newpupil;

% zscore since we work with the bandpassed signal
data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:) = ...
    zscore(data.trial{1}(find(strcmp(data.label, 'EyePupil')==1),:));

% ==================================================================
% define trials
% ==================================================================

cfg                         = [];
cfg.dataset                 = ascFile;
cfg.event                   = event;
cfg.trialfun                = 'my_trialfun';
cfg.trialdef.pre            = 0;
cfg.trialdef.post           = 2;
cfg.fsample                 = asc.fsample;
cfg.sj                      = 22;
cfg.session                 = 5;
cfg.block                   = 1;
cfg                         = ft_definetrial(cfg);
% epoch
data                        = ft_redefinetrial(cfg, data);  

% note that this trialfun defines four events per trial: two stimuli (that
% have to be compared), a response, and feedback). We have the trial start
% at fixation and save the samples of each of those events (as well as more
% information about stimulus identity, choice and accuracy) in
% data.trialinfo:
% [fixoffset refoffset stimtype stimoffset resptype respcorrect respoffset ...
% feedbacktype feedbackoffset trlcnt blockcnt session]
  
% ==================================================================
% downsample before saving
% ==================================================================

cfg             = [];
cfg.resamplefs  = 100;
cfg.fsample     = data.fsample;

% see Niels' message on the FT mailing list
samplerows = find(data.trialinfo(1,:)>100); % indices of the rows with sample values (and not event codes)
data.trialinfo(:,samplerows) = round(data.trialinfo(:,samplerows) * (cfg.resamplefs/cfg.fsample));

% use fieldtrip to resample
data    = ft_resampledata(cfg, data);

% ==================================================================
% visualize the pupil timecourse
% ==================================================================

cfg                 = [];
cfg.channel         = 'EyePupil';
cfg.trials(1).name  = 'correct';
cfg.trials(1).idx   = find(data.trialinfo(:, 8) == 1);
cfg.trials(2).name  = 'error';
cfg.trials(2).idx   = find(data.trialinfo(:, 8) == 0);

clf; subplot(2,2,[1 2]);
plotEventRelated(cfg, data);
title(ascFile, 'interpreter', 'none');
ylabel('Pupil response (z)');

% ==================================================================
% save file
% ==================================================================

filename = regexprep(ascFile, 'asc', 'mat');
save(filename, 'data');
