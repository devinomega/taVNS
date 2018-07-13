function varargout = respPick(varargin)
% RESPPICK MATLAB code for respPick.fig
%      RESPPICK, by itself, creates a new RESPPICK or raises the existing
%      singleton*.
%
%      H = RESPPICK returns the handle to a new RESPPICK or the handle to
%      the existing singleton*.
%
%      RESPPICK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESPPICK.M with the given input arguments.
%
%      RESPPICK('Property','Value',...) creates a new RESPPICK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before respPick_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to respPick_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help respPick

% Last Modified by GUIDE v2.5 15-Jun-2018 18:38:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @respPick_OpeningFcn, ...
                   'gui_OutputFcn',  @respPick_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before respPick is made visible.
function respPick_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to respPick (see VARARGIN)

handles.break = 0;  %if someone closes the figure then you want to break out of ECG_preProc
set(0,'userdata',{[] []});   %reset the root userdata

% Check for correct inputs
if nargin< 4    
    plot(rand(5));  %plot something random if not input
else 
    handles.rawData = varargin{1};
    handles.rawTime  = (0:numel(handles.rawData)-1);
    handles.secTime  = (0:numel(handles.rawData)-1)/2000;
    
    baseLnCor = detrend(handles.rawData);
    %filter data
    highEnd = 1;
    f_s = 2000; %sampling frequency
    filterOrder = 2; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
%     [b, a] = butter(filterOrder, [0.5 highEnd]/(f_s/2),'low'); % Generate filter coefficients
    [b, a] = butter(filterOrder, [0.05 highEnd]/(f_s/2)); % Generate filter coefficients
    handles.filtData = filtfilt(b, a, baseLnCor); % Apply filter to data using zero-phase filtering
    
    handles.threshold =  0;
    [handles.peaksMax,handles.locsMax] = findpeaks(handles.filtData,handles.rawTime,'MinPeakDistance',3000);
    [handles.peaksMin,handles.locsMin] = findpeaks(handles.filtData*-1,handles.rawTime,'MinPeakDistance',3000);
    handles.peaksMin = handles.peaksMin*-1;
    
    % Plot the filtered detrended data
    handles.ax(1) = handles.transform_graph;
    plot(handles.transform_graph,handles.secTime, handles.filtData,'Parent',handles.ax(1))
    title(handles.transform_graph,'(In/Ex)halation Localized by findPeaks')
    dcm_obj = datacursormode(gcf);
    datacursormode on
    set(dcm_obj,'UpdateFcn',@myupdatefcn);
    hold(handles.transform_graph,'on')
    
%   Plot the inhilation 
    plot(handles.transform_graph,handles.secTime(handles.locsMax),handles.peaksMax,'ro','Parent',handles.ax(1))
    %Plot the exhalation
    plot(handles.transform_graph,handles.secTime(handles.locsMin),handles.peaksMin,'mo','Parent',handles.ax(1))
    xlabel(handles.transform_graph,'Seconds')
    hold(handles.transform_graph,'off')
%     
    %Plot the raw respiration signal   
    handles.ax(2) = handles.rawEEG_graph;
    plot(handles.rawEEG_graph, handles.secTime,handles.rawData,'Parent',handles.ax(2))
    title(handles.rawEEG_graph,'Raw resp')
    xlabel(handles.rawEEG_graph,'Seconds')
    dcm_obj = datacursormode(gcf);  %may not work
    datacursormode on
    set(dcm_obj,'UpdateFcn',@myupdatefcn);
    hold(handles.rawEEG_graph,'on')
    
    %   Plot the inhilation 
    plot(handles.rawEEG_graph,handles.secTime(handles.locsMax),handles.rawData(handles.locsMax),'ro','Parent',handles.ax(2))
    %Plot the exhalation
    plot(handles.rawEEG_graph,handles.secTime(handles.locsMin),handles.rawData(handles.locsMin),'mo','Parent',handles.ax(2))
    xlabel(handles.rawEEG_graph,'Seconds')
    hold(handles.rawEEG_graph,'off')
    
    %Link the axes for zoom and turn on the tool bar
    linkaxes(handles.ax,'x')
    set(gcf,'toolbar','figure');    %Enable zoom and data cursor (all I care about)
    
    handles.AddPoints = [];   %The points added manually
%     pos = get(0,'userdata');
end

%If there's a subject name - put it in the title
if numel(varargin) == 2
    set(handles.subjName,'String',varargin{2})
end

% Choose default command line output for respPick
handles.output = {handles.locsMin handles.locsMax};
set(handles.threshBox,'String',num2str(handles.threshold))
set(handles.PointTable,'Data',[sort([handles.secTime(handles.locsMin) handles.secTime(handles.locsMax)])' [60./diff(sort([handles.secTime(handles.locsMin) handles.secTime(handles.locsMax)])) NaN]'])

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RR_Peak_Pick wait for user response (see UIRESUME)
uiwait(handles.figure1);
% UIWAIT makes respPick wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = respPick_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.break
    handles.output = {999 999};
end
% Get default command line output from handles structure
varargout{1} = handles.output{1};
varargout{2} = handles.output{2};
% The figure can be deleted now
delete(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%% CloseReqFCN
%turn on later
if isequal(get(hObject, 'waitstatus'), 'waiting')
    handles.break = 1;
    % The GUI is still in UIWAIT, use UIRESUME
    uiresume(hObject);
else
    handles.break = 1;
    % The GUI is no longer waiting, just close it
    delete(hObject);
end
guidata(hObject,handles)

% --- Executes on button press in AddPnt.
function AddPnt_Callback(hObject, eventdata, handles)
% hObject    handle to AddPnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.transform_graph);
xlim = get(handles.ax(1), 'XLim');
ylim = get(handles.ax(1),'YLim');
ylim_ECG = get(handles.ax(2),'YLim');
cla;

pntData = get(0,'userdata');  %Get the input from the graph
cPos = pntData{1};   %all x positions
cFig = pntData{2};  %which figure it came from 
handles.AddPoints =[handles.AddPoints;cPos];
curData = get(handles.PointTable,'Data');   %get all the current pointtable data

for p = 1:numel(cFig)
    tempFig = cFig{p};
    pos = cPos(p,:);
    % determine which graph the point was drawn from
    if strcmp(tempFig,'Raw ECG')
        %If it was the raw - then the peak needs to be identified on the
        %filtered
        addPeak= handles.filtData(handles.secTime==pos(1));
    else
        %otherwise keep the value
        addPeak = pos(2);
    end

    newData = sort([curData(:,1); pos(1)]); %add the new position and sort it
    curData = [newData [(60./diff(newData));NaN]];  %recalculate the BPM    
    nn = curData(logical(abs(diff(curData(:,1)<pos(1)))),1);   %identify the point just prior to the one added
    
    % If the point prior was a Max - add new point as a Min, and vice versa
    if sum(handles.secTime(handles.locsMax) == nn)
        [handles.locsMin, Indx] = sort([handles.locsMin find(handles.secTime==pos(1))]);
        handles.peaksMin = [handles.peaksMin addPeak];    %sort the peaks as well
        handles.peaksMin = handles.peaksMin(Indx);
    else
        [handles.locsMax, Indx] = sort([handles.locsMax find(handles.secTime==pos(1))]);
        handles.peaksMax = [handles.peaksMax addPeak];    %sort the peaks as well
        handles.peaksMax = handles.peaksMax(Indx);
    end
    
end

% Replot the filtered data
plot(handles.transform_graph,handles.secTime, handles.filtData,'Parent',handles.ax(1))
zoom reset
set(handles.ax(1),'XLim',xlim)
set(handles.ax(1),'YLim',ylim)
title('(In/Ex)halation  Localized by findPeaks')
hold(handles.transform_graph,'on')

% Replot peak points
plot(handles.transform_graph,handles.secTime(handles.locsMax),handles.peaksMax,'ro')
plot(handles.transform_graph,handles.secTime(handles.locsMin),handles.peaksMin,'mo')
xlabel('Seconds')
hold off

axes(handles.rawEEG_graph);

%Plot the raw Resp signal   
plot(handles.rawEEG_graph, handles.secTime,handles.rawData,'Parent',handles.ax(2))
zoom reset
set(handles.ax(2),'XLim',xlim)
set(handles.ax(2),'YLim',ylim_ECG)
title(handles.rawEEG_graph,'Raw Resp')
hold(handles.rawEEG_graph,'on')

% Plot the Peaks for the raw Resp data
plot(handles.rawEEG_graph,handles.secTime(handles.locsMax),handles.rawData(handles.locsMax),'ro','Parent',handles.ax(2))
plot(handles.rawEEG_graph,handles.secTime(handles.locsMin),handles.rawData(handles.locsMin),'mo','Parent',handles.ax(2))
xlabel(handles.rawEEG_graph,'Seconds')
hold(handles.rawEEG_graph,'off')

set(handles.PointTable,'Data',curData)  %Set the point table again
set(0,'userdata',{[] []});  %reset the user data 
handles.output = {handles.locsMin handles.locsMax}; %Update the output

guidata(hObject, handles);

% --- Executes on button press in Update.
function Update_Callback(hObject, eventdata, handles)
% hObject    handle to Update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.transform_graph);
cla;

%Update the trehsold and replot it
handles.threshold = get(handles.threshBox,'String');    %turn back on
[handles.peaksMax,handles.locsMax] = findpeaks(handles.filtData,handles.rawTime,'MinPeakDistance',3000);
[handles.peaksMin,handles.locsMin] = findpeaks(handles.filtData*-1,handles.rawTime,'MinPeakDistance',3000);
handles.peaksMin = handles.peaksMin*-1;

%Add back the points that were manually entered
if ~isempty(handles.AddPoints)
    if pos(1) > handles.threshold
        [handles.locsMax, Indx] = sort([handles.locsMax handles.AddPoints(1)]);
        handles.peaksMax = [handles.peaksMax handles.AddPoints(2)];
        handles.peaksMax = handles.peaksMax(Indx);
    else
        [handles.locsMin, Indx] = sort([handles.locsMin handles.AddPoints(1)]);
        handles.peaksMin = [handles.peaksMin handles.AddPoints(2)];
        handles.peaksMin = handles.peaksMin(Indx);
    end
end

%Plot the filtered data
plot(handles.transform_graph,handles.secTime, handles.filtData)
title('(In/Ex)halation  Localized by findPeaks')
hold(handles.transform_graph,'on')

% Replot peak points
plot(handles.transform_graph,handles.secTime(handles.locsMax),handles.peaksMax,'ro')
plot(handles.transform_graph,handles.secTime(handles.locsMin),handles.peaksMin,'mo')
xlabel('Seconds')
hold off

axes(handles.rawEEG_graph);

%Plot the raw Respsignal   
plot(handles.rawEEG_graph, handles.secTime,handles.rawData,'Parent',handles.ax(2))
title(handles.rawEEG_graph,'Raw resp')
xlabel(handles.rawEEG_graph,'Seconds')
dcm_obj = datacursormode(gcf);  %may not work
datacursormode on
set(dcm_obj,'UpdateFcn',@myupdatefcn);
hold(handles.rawEEG_graph,'on')

%   Plot the Peaks for the raw ECG data
plot(handles.rawEEG_graph,handles.secTime(handles.locsMax),handles.rawData(handles.locsMax),'ro','Parent',handles.ax(2))
plot(handles.rawEEG_graph,handles.secTime(handles.locsMin),handles.rawData(handles.locsMin),'mo','Parent',handles.ax(2))
xlabel(handles.rawEEG_graph,'Seconds')
hold(handles.rawEEG_graph,'off')

%redraw the point table
set(handles.PointTable,'Data',[sort([handles.locsMin handles.locsMax])' [60./diff(sort([handles.locsMin handles.locsMax])) NaN]'])

%Update output
handles.output = {handles.locsMin handles.locsMax};

guidata(hObject, handles);


function threshBox_Callback(hObject, eventdata, handles)
% hObject    handle to threshBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshBox as text
%        str2double(get(hObject,'String')) returns contents of threshBox as a double

get(hObject,'String'); %returns contents of threshBox as text
str2double(get(hObject,'String')); %returns contents of threshBox as a double
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function threshBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Accept.
function Accept_Callback(hObject, eventdata, handles)
% hObject    handle to Accept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles)
uiresume(handles.figure1)

% --- Executes on button press in Reject.
function Reject_Callback(hObject, eventdata, handles)
% hObject    handle to Reject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.locs = NaN;

%Update output
handles.output = handles.locs;
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in RejectPts.
function RejectPts_Callback(hObject, eventdata, handles)
% hObject    handle to RejectPts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.transform_graph);
xlim = get(handles.ax(1), 'XLim');
ylim = get(handles.ax(1),'YLim');
ylim_ECG = get(handles.ax(2),'YLim');
cla;
eFlag = 1;

curData = get(handles.PointTable,'Data');
pntData = get(0,'userdata');
%Get the rows that need to be deleted
try
    rows = handles.rejectPtsData.Indices;
catch ME
    if (strcmp(ME.identifier,'MATLAB:nonExistentField'))
    	if ~isempty(pntData{1})
            rows = pntData{1}(:,1);   %all x positions
            rows = dsearchn(curData(:,1),rows);
%         cFig = pntData{2};  %which figure it came from 
        else
            eFlag = 0;
        end
    else
        eFlag = 0;  %no points selected and 
    end
end

if eFlag
    indxMax = logical(sum(handles.secTime(handles.locsMax) == curData(rows(:,1),1),1));
    indxMin = logical(sum(handles.secTime(handles.locsMin) == curData(rows(:,1),1),1));

    if sum(indxMax) >= 1
    % Update figure and table
    handles.locsMax(indxMax) = [];   %Locs w/ point removed
    handles.peaksMax(indxMax) = [];  %amp of point removed
    end
    if sum(indxMin) >= 1
        % Update figure and table
        handles.locsMin(indxMin) = [];   %Locs w/ point removed
        handles.peaksMin(indxMin) = [];  %amp of point removed
    end
    % BPMNew = [60./diff([handles.locsMin handlesMax]) NaN'];%BPM recalculated
    % curData = [[handles.locsMin locsMax]' BPMNew];   %
    curData(rows(:,1),:) = [];
    set(handles.PointTable,'Data',curData)
    
    % Plot filter Data
    plot(handles.transform_graph,handles.secTime, handles.filtData)
    title('(In/Ex)halation  Localized by findPeaks')
    hold(handles.transform_graph,'on')
    
    %replot peak points
    plot(handles.transform_graph,handles.secTime(handles.locsMax),handles.peaksMax,'ro')
    plot(handles.transform_graph,handles.secTime(handles.locsMin),handles.peaksMin,'mo')
    % zoom on
    xlabel('Seconds')
    hold off
    
    axes(handles.rawEEG_graph);
    %Plot the raw Respsignal
    plot(handles.rawEEG_graph, handles.secTime,handles.rawData,'Parent',handles.ax(2))
    title(handles.rawEEG_graph,'Raw resp')
    xlabel(handles.rawEEG_graph,'Seconds')
    dcm_obj = datacursormode(gcf);  %may not work
    datacursormode on
    set(dcm_obj,'UpdateFcn',@myupdatefcn);
    hold(handles.rawEEG_graph,'on')
    
    %   Plot the Peaks for the raw ECG data
    plot(handles.rawEEG_graph,handles.secTime(handles.locsMax),handles.rawData(handles.locsMax),'ro','Parent',handles.ax(2))
    plot(handles.rawEEG_graph,handles.secTime(handles.locsMin),handles.rawData(handles.locsMin),'mo','Parent',handles.ax(2))
    xlabel(handles.rawEEG_graph,'Seconds')
    hold(handles.rawEEG_graph,'off')
    
    %Update output
    handles.output = {handles.locsMin handles.locsMax};
    set(0,'userdata',{[] []});  %reset the point selection
end
guidata(hObject,handles)

% --- Executes when selected cell(s) is changed in PointTable.
function PointTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to PointTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
handles.rejectPtsData = eventdata;
guidata(hObject, handles);

function [txt, curFig] = myupdatefcn(~, event_obj)
    temp = get(0,'userdata');
    %Fix - if point has previously been selected, remove from list
	unslctPnt = (event_obj.Position == temp{1});
    if sum(event_obj.position == temp{1})
        temp{1}(unslctPnt) = [];
        temp{2}(unslctPnt) = [];
        pos = temp{1};
        curFig = temp{2};
    else
        pos = [temp{1};event_obj.Position];
        a = get(gca,'Title');
    %   disp(['You clicked X: ',num2str(handles.pos(1)),', Y: ',num2str(handles.pos(2))]);
        txt = pos(1);
        curFig = {get(a,'String');temp{2}};
    end
        set(0,'userdata',{pos, curFig});
