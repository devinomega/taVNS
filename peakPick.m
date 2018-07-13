function varargout = peakPick(varargin)
% PEAKPICK MATLAB code for peakPick.fig
%      PEAKPICK, by itself, creates a new PEAKPICK or raises the existing
%      singleton*.
%
%      H = PEAKPICK returns the handle to a new PEAKPICK or the handle to
%      the existing singleton*.
%
%      PEAKPICK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PEAKPICK.M with the given input arguments.
%
%      PEAKPICK('Property','Value',...) creates a new PEAKPICK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before peakPick_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to peakPick_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help peakPick

% Last Modified by GUIDE v2.5 03-May-2018 10:41:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @peakPick_OpeningFcn, ...
                   'gui_OutputFcn',  @peakPick_OutputFcn, ...
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


% --- Executes just before peakPick is made visible.
function peakPick_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to peakPick (see VARARGIN)

handles.break = 0;  %if someone closes the figure then you want to break out of ECG_preProc

% Check for correct inputs
if nargin< 4    
    plot(rand(5));  %plot something random if not input
else 
    handles.rawData = varargin{1};
    handles.rawTime = 0:numel(handles.rawData)-1;
    handles.secTime  = (0:numel(handles.rawData)-1)/2000;
    
    %filter data
    lowEnd = .5;
    highEnd = 26;
    f_s = 2000; %sampling frequency
    
    filterOrder = 2; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
    [b, a] = butter(filterOrder, [lowEnd highEnd]/(f_s/2)); % Generate filter coefficients
    handles.filtData = filtfilt(b, a,handles.rawData); % Apply filter to data using zero-phase filtering
    
    % use wavelet isolate the R peaks 
    %   Use findpeaks to determine the peak locations. 
    %   Plot the R-peak waveform along with the expert and automatic annotations
    % db6
%     [e,f]=wavedec(handles.filtData,10,'db6');% Wavelet implementation
%     g=wrcoef('a',e,f,'db6',8); 
%     handles.wavelet=handles.filtData-g; % subtracting 10th level aproximation signal
%                        %from original signal 
%     handles.wavelet(handles.wavelet<0) = 0;
%     handles.wavelet = handles.wavelet.^2;
%     handles.threshold =  quantile(handles.wavelet,.98);

    % sym8 - no filtering
    wt = modwt(handles.rawData,'sym8');
    wtrec = zeros(size(wt));
    wtrec(6:7,:) = wt(6:7,:);
    handles.wavelet = imodwt(wtrec,'sym8');
    handles.wavelet(handles.wavelet<-0) = 0;
    handles.wavelet = handles.wavelet.^2;
    handles.threshold =  quantile(handles.wavelet,.98);
    
    % sym4 - no filtering
%     wt = modwt(handles.rawData,5);
%     wtrec = zeros(size(wt));
%     wtrec(4:5,:) = wt(4:5,:);
%     handles.wavelet = imodwt(wtrec,'sym4');
%     handles.wavelet(handles.wavelet<-0) = 0;
%     handles.wavelet = handles.wavelet.^2;
%     handles.threshold =  quantile(handles.wavelet,.98);

    [handles.peaks,handles.locs] = findpeaks(handles.wavelet,handles.rawTime,'MinPeakHeight',handles.threshold,'MinPeakDistance',1000);
    
    % Plot the Wavelet data
    handles.ax(1) = handles.transform_graph;
    plot(handles.transform_graph,handles.secTime, handles.wavelet,'Parent',handles.ax(1));
    title(handles.transform_graph,'R-Waves Localized by Wavelet Transform')
    dcm_obj = datacursormode(gcf);
    datacursormode on
    set(dcm_obj,'UpdateFcn',@myupdatefcn);
    hold(handles.transform_graph,'on')
    
%   Plot the Peaks for the wavelet data
    plot(handles.transform_graph,handles.secTime(handles.locs),handles.peaks,'ro','Parent',handles.ax(1))
    xlabel(handles.transform_graph,'Seconds')
    hold(handles.transform_graph,'off')
    
    %Plot the raw ECG signal   
    handles.ax(2) = handles.rawEEG_graph;
    plot(handles.rawEEG_graph, handles.secTime,handles.rawData,'Parent',handles.ax(2))
    title(handles.rawEEG_graph,'Raw ECG')
    xlabel(handles.rawEEG_graph,'Seconds')
    hold(handles.rawEEG_graph,'on')
    
%   Plot the Peaks for the raw ECG data
    plot(handles.rawEEG_graph,handles.secTime(handles.locs),handles.rawData(handles.locs),'ro','Parent',handles.ax(2))
    xlabel(handles.rawEEG_graph,'Seconds')
    hold(handles.rawEEG_graph,'off')
    
    %Link the axes for zoom and turn on the tool bar
    linkaxes(handles.ax,'x')
    set(gcf,'toolbar','figure');    %Enable zoom and data cursor (all I care about)
    
    handles.AddPoints = [];   %The points added manually
end

%If there's a subject name - put it in the title
if numel(varargin) == 2
    set(handles.subjName,'String',varargin{2})
end

% Choose default command line output for RR_Peak_Pick
handles.output = handles.locs;
set(handles.threshBox,'String',num2str(handles.threshold))
set(handles.PointTable,'Data',[handles.secTime(handles.locs)' [60./diff(handles.secTime(handles.locs)) NaN]'])

% set(gcf, 'units', 'normalized', 'position', [0 0 1 1])

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RR_Peak_Pick wait for user response (see UIRESUME)
uiwait(handles.figure1);
% UIWAIT makes peakPick wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = peakPick_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.break
    handles.output = 999;
end
% Get default command line output from handles structure
varargout{1} = handles.output;
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

pntData = get(0,'userdata');  %turn this on 
pos = pntData{1};
cFig = pntData{2};
handles.AddPoints =[handles.AddPoints;pos];

if strcmp(cFig,'Raw ECG')
    addPeak= handles.wavelet(handles.secTime==pos(1));
else
    addPeak = pos(2);
end

% Check this
curData = get(handles.PointTable,'Data');
newData = sort([curData(:,1); pos(1)]);
curData = [newData [(60./diff(newData));NaN]];
set(handles.PointTable,'Data',curData)

% Update figure as well
[handles.locs, Indx] = sort([handles.locs find(handles.secTime==pos(1))]);
handles.peaks = [handles.peaks addPeak];
handles.peaks = handles.peaks(Indx);    %why???

%replot the figure
plot(handles.transform_graph,handles.secTime, handles.wavelet,'Parent',handles.ax(1))
zoom reset
set(handles.ax(1),'XLim',xlim)
set(handles.ax(1),'YLim',ylim)
title('R-Waves Localized by Wavelet Transform')
hold(handles.transform_graph,'on')

%replot peak points
plot(handles.transform_graph,handles.secTime(handles.locs),handles.peaks,'ro')
xlabel('Seconds')
hold off

axes(handles.rawEEG_graph);
%Plot the raw ECG signal   
plot(handles.rawEEG_graph, handles.secTime,handles.rawData,'Parent',handles.ax(2))
zoom reset
set(handles.ax(2),'XLim',xlim)
set(handles.ax(2),'YLim',ylim_ECG)
title(handles.rawEEG_graph,'Raw ECG')
hold(handles.rawEEG_graph,'on')

%   Plot the Peaks for the raw ECG data
plot(handles.rawEEG_graph,handles.secTime(handles.locs),handles.rawData(handles.locs),'ro','Parent',handles.ax(2))
xlabel(handles.rawEEG_graph,'Seconds')
hold(handles.rawEEG_graph,'off')

%Update the output
handles.output = handles.locs;

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
[handles.peaks,handles.locs] = findpeaks(handles.wavelet,handles.rawTime,'MinPeakHeight',str2double(handles.threshold),'MinPeakDistance',1000);

%Add back the points that were manually entered
if ~isempty(handles.AddPoints)
    [handles.locs, Indx] = sort([handles.locs handles.rawTime(handles.secTime==handles.AddPoints(:,1))]);
    handles.peaks = [handles.peaks handles.AddPoints(:,2)];
    handles.peaks = handles.peaks(Indx);
end

plot(handles.transform_graph,handles.secTime, handles.wavelet)
% zoom on
title('R-Waves Localized by Wavelet Transform')
hold(handles.transform_graph,'on')

plot(handles.transform_graph,handles.secTime(handles.locs),handles.peaks,'ro')
% zoom on
xlabel('Seconds')
hold off

axes(handles.rawEEG_graph);
%Plot the raw ECG signal   
plot(handles.rawEEG_graph, handles.secTime,handles.rawData,'Parent',handles.ax(2))
title(handles.rawEEG_graph,'Raw ECG')
xlabel(handles.rawEEG_graph,'Seconds')
dcm_obj = datacursormode(gcf);  %may not work
datacursormode on
set(dcm_obj,'UpdateFcn',@myupdatefcn);
hold(handles.rawEEG_graph,'on')

%   Plot the Peaks for the raw ECG data
plot(handles.rawEEG_graph,handles.secTime(handles.locs),handles.rawData(handles.locs),'ro','Parent',handles.ax(2))
xlabel(handles.rawEEG_graph,'Seconds')
hold(handles.rawEEG_graph,'off')

set(handles.PointTable,'Data',[handles.secTime(handles.locs)' [60./diff(handles.secTime(handles.locs)) NaN]'])

%Update output
handles.output = handles.locs;

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

%Get the rows that need to be deleted
rows = handles.rejectPtsData.Indices;

% Update figure and table
handles.locs(rows(:,1)) = [];   %Locs w/ point removed
handles.peaks(rows(:,1)) = [];  %amp of point removed
BPMNew = [60./diff(handles.secTime(handles.locs)) NaN]';  %BPM recalculated
curData = [handles.secTime(handles.locs)' BPMNew];   %
set(handles.PointTable,'Data',curData)

% Plot graphs
plot(handles.transform_graph,handles.secTime, handles.wavelet,'Parent',handles.ax(1))
zoom reset
set(handles.ax(1),'XLim',xlim)
set(handles.ax(1),'YLim',ylim)
title('R-Waves Localized by Wavelet Transform')
hold(handles.transform_graph,'on')

plot(handles.transform_graph,handles.secTime(handles.locs),handles.peaks,'ro')
xlabel('Seconds')
hold off

axes(handles.rawEEG_graph);
%Plot the raw ECG signal   
plot(handles.rawEEG_graph, handles.secTime,handles.rawData,'Parent',handles.ax(2))
zoom reset
set(handles.ax(2),'XLim',xlim)
set(handles.ax(2),'YLim',ylim_ECG)

title(handles.rawEEG_graph,'Raw ECG')
hold(handles.rawEEG_graph,'on')

%   Plot the Peaks for the raw ECG data
plot(handles.rawEEG_graph,handles.secTime(handles.locs),handles.rawData(handles.locs),'ro','Parent',handles.ax(2))
xlabel(handles.rawEEG_graph,'Seconds')
hold(handles.rawEEG_graph,'off')

%Update output
handles.output = handles.locs;
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
    pos = event_obj.Position;
    a = get(gca,'Title');
    %   disp(['You clicked X: ',num2str(handles.pos(1)),', Y: ',num2str(handles.pos(2))]);
    txt = pos(1);
    curFig = get(a,'String');
    set(0,'userdata',{pos, curFig});
