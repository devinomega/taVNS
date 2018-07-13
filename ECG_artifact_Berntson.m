%% Start Program
% 1) Check for criterion level beat differences
% 2) If a criterion difference is found, an artifact flag is set and the
%   integrity of the prior beat difference is evaluated
% 3) If the prio beat difference is within criterion, the program then
%   determines whether the criterion-level difference was due to the presence
%   of a long beat or a short beat and control is then passed to the
%   appropriate subroutine
%
% Criterion Beat Difference = (MAD+MED)/2
%
% MAD = Minimum Artificat Deviation
%       Shortest expected veridical beat (SEB)/3
%       SEB = Median Beat - 2.9*Quartile Deviation
%       Quartile Deviation = Q3-Q1/2
% MED = Maximum Expected Deviation of veridical beats (Like STD)
%       MED = +/- 3.32 QD
% MAD> MED index *ideally

defaultPath = 'C:\Users\devinomega\Dropbox\MATLAB\tVNS\NI\subjData\';
subj = 'BIGMAT_TV_20_02.mat';
load([defaultPath subj])
temp = ECG.locsData{2};
% temp = newLocs;
rawData = ECG.allData{2}(1:(2000*135));

diffLocs =diff(temp*1000);
tm  = (0:(2000*135)-1)/2000;

%% Criterion
% MED
QD = (quantile(diffLocs,0.75)-quantile(diffLocs,0.25))/2;   %Quartile deviaiton
MED = 3.32*QD;  %Minimum expected difference
% meanDiff = diffLocs-mean(diffLocs);    %QD score all the data 

% MAD
SEB = median(diffLocs)-(QD*2.9);    %shortest expected veridical beat
MAD = SEB/3; % Minimum artifact difference

CBD = (MED+MAD)/2;  %Criterion beat difference

% outliers = abs(meanDiff)>CBD;
beatDiff = diff(diffLocs);
outliers = abs(beatDiff)>CBD;
%%
%Flag arrays
artifact = outliers;    %correctly identified artifact
shrtLng =  zeros(size(outliers));    %whether the beat was long (1) or short (2)
FA = zeros(size(outliers)); %false alarm

indx = 1:numel(outliers);   %index of the outliers
indx = indx(outliers);

for n = 1:sum(outliers)
    % Is it beat (RR Diff) to beat (RR diff) differences that are set to
    % the criterion? Or is beat with the mean subtracted compare to the
    % criterion?
    curIndx = indx(n);

        
        if abs(beatDiff(curIndx-1)) < CBD   %If the previous beat is not an artifact
            %check to see if the current artifact is a short beat (>0) or long (<0)
            if (beatDiff(curIndx)*-1)  < 0  %B(n)-B(n+1)
                %case for long beat
                if abs(beatDiff(curIndx+2)) < CBD % are the beats two ahead not artifacts?
                    X = diffLocs(curIndx+1)/2;
                    if ((X-diffLocs(curIndx)) < (-1*CBD)) && ((X-diffLocs(curIndx+2)) < (-1*CBD))
                        FA(curIndx) = 1;
                        artifact(curIndx) = 0;
                    else
                        shrtLng(curIndx) = 1;
                    end
                else
                   shrtLng(curIndx) = 3;
                end
            else
                if abs(beatDiff(curIndx+2))< CBD % are the beat two ahead not artifacts?
                % case for short beat
                    if diffLocs(curIndx)< diffLocs(curIndx+2)
                        X = diffLocs(curIndx) + diffLocs(curIndx+1);
                    else
                        X = diffLocs(curIndx+1) + diffLocs(curIndx+2); 
                    end
                
                    if X-diffLocs(curIndx-1) > CBD && X-diffLocs(curIndx+1) > CBD 
                        FA(curIndx) = 1;
                        artifact(curIndx) = 0;
                    else
                        shrtLng(curIndx) = 2;
                    end
                else
                    shrtLng(curIndx) = 3;    
                end
            end
        else  %If the previous beat was an artifact
            shrtLng(curIndx) = 3;   %can't evaluate
        end
end
%do the calculations to either add a point or remove a point (see which one
%it is by adding it to different arrays)
%is flagged (don't do it again on FA)
%Rerun the criterion check - if there are still artifacts
%do some stuff down here where either flag is set and you're all done - o

%                        newLocs = [newLocs(1:(curIndx+1)) (newLocs(curIndx+1)+(X/1000)) newLocs((curIndx+2):end)];
% f = 1; %while

artLoc = find(artifact);
%Check to make sure everything is copasetic
figure
plot(tm,rawData)
hold on
% plot(temp,rawData(ismember(tm,temp)),'or')
plot(temp,rawData(dsearchn(tm',temp')),'or');
% newLocs = temp;
clrs = {[0.6 0.4 0.9] [0 1 0] [1 0 0]};
for p = 1:sum(artifact)
    ptchidx= [temp(artLoc(p)+1)  temp(artLoc(p)+2)];
    patch([ptchidx fliplr(ptchidx)],reshape(repmat(ylim,[2 1]),[1 4]),clrs{shrtLng(artLoc(p))}, 'FaceAlpha',0.3, 'EdgeColor','none')
end

newLocs = temp;
% newLocs = [newLocs(1:(curIndx+1)) (newLocs(curIndx+1)+(X/1000)) newLocs((curIndx+2):end)];
lngPts =newLocs(find((shrtLng == 1)) + 1) + (diffLocs(find((shrtLng == 1)) + 1)./2)/1000;
newLocs = sort([newLocs lngPts]);