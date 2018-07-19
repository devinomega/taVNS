sampleLen = 135*2000;
sampleRate = 2000;
%%
Fpath = ('C:\Users\bikson\Documents\taVNS_data\testSubjData\');
FileHolder = dir('C:\Users\bikson\Documents\taVNS_data\testSubjData\BIGMAT*.mat');     
for i = 1:numel(FileHolder)    
  curSubj = FileHolder(i).name; %subject currently being worked on
  comboFile =  [Fpath curSubj]; %combines to make one string 
  load(comboFile)
                       
%%
    B = zeros(5,3);
for  block = 1:5 
    Beep = RespRate.locsMin{block};
    if sum(mod(Beep,1) ~= 0) > 0
        Beep = Beep*2000;
    end
    MaxDiff = diff(Beep);
    MaxSec = 1./(MaxDiff/2000); % setting data to resp per min
    Max15 = (Beep/2000)<15;     % first set of 15 points
    Avg1 = mean(Max15);
    Max60int = ((Beep/2000)>15) & ((Beep/2000)<75);  % middle set between 15 and 75
    Avg2 = mean(Max60int);
    Max60Final = (Beep/2000)>75;  % final set of numbers after 75
    Avg3 = mean(Max60Final);
    
    endData = [Avg1 Avg2 Avg3];
    B(block,:) = endData;
    
end
RespRate.AveragesMin = B;
save(comboFile,'RespRate','-append')
end