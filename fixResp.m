% aPath = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\072018_subjData\'; %postData saved in locs
% bPath=  'C:\Users\devinomega\Dropbox\MATLAB\tVNS\071118_subjData\';   %correct locs data
% sPath=  'C:\Users\devinomega\Dropbox\MATLAB\tVNS\subjData\';   %save file location
% 
% aFiles = dir([aPath '\BIGMAT*.mat']);
% 
% for i = 1:numel(aFiles)
%     curSubj = aFiles(i).name; %subject currently being worked on
%     aData = load([aPath curSubj]);   %postData saved in locs
%     bData = load([bPath curSubj],'RespRate');   %correct locs data, no postLocs
%     
%     aData.RespRate.postMin = aData.RespRate.locsMin;
%     aData.RespRate.postMax = aData.RespRate.locsMax;
%     
%     aData.RespRate.locsMin = bData.RespRate.locsMin;
%     aData.RespRate.locsMax = bData.RespRate.locsMax;
%     
%     fields = {'min','max'};
%     RespRate = rmfield(aData.RespRate,fields);
%     
%     copyfile([aPath curSubj],[sPath curSubj])
%     save([sPath curSubj],'RespRate','-append')
% end

aPath = 'C:\Users\devinomega\Desktop\subjOnlyPostBase\'; %postData saved in locs
% bPath=  'C:\Users\devinomega\Dropbox\MATLAB\tVNS\071118_subjData\';   %correct locs data
sPath=  'C:\Users\devinomega\Dropbox\MATLAB\tVNS\subjData\';   %save file location

aFiles = dir([aPath '\BIGMAT*.mat']);

for i = 1:numel(aFiles)
    curSubj = aFiles(i).name; %subject currently being worked on
    aData = load([aPath curSubj]);   %postData saved in locs
%     bData = load([bPath curSubj],'RespRate');   %correct locs data, no postLocs
    
    aData.RespRate.postMin = aData.RespRate.locsMin;
    aData.RespRate.postMax = aData.RespRate.locsMax;
    
    aData.RespRate.locsMin = [];
    aData.RespRate.locsMax = [];
    
    fields = {'min','max'};
    RespRate = rmfield(aData.RespRate,fields);
    
    copyfile([aPath curSubj],[sPath curSubj])
    save([sPath curSubj],'RespRate','-append')
end