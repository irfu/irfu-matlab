function varargout = tt_browse(varargin)
% IRF.TT_BROWSE browse time table
%	IRF.TT_BROWSE(TT) browse time table TT
%		time table shoud be in irf.TimeTable format
%
% 	IRF.TT_BROWSE(TT,'plotfile') browse and make figure defined in plotfile
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% $Revision$  $Date$

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       '+irf/tt_browse.fig', ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @tt_browse_OpeningFcn, ...
  'gui_OutputFcn',  @tt_browse_OutputFcn, ...
  'gui_LayoutFcn',  [] , ...
  'gui_Callback',   []);
if nargin && ischar(varargin{1})
  gui_State.gui_Callback = str2func(varargin{1});
elseif nargin == 0
  help irf.tt_browse;
  return;
end

if nargout
  [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
  gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before irf.tt_browse is made visible.
function tt_browse_OpeningFcn(hObject,~,handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to irf.tt_browse (see VARARGIN)

% Choose default command line output for IRF.TT_BROWSE
handles.output = hObject;
handles.workDirectory=pwd;
isSpecifiedDataDirectory = 0;
if isa(varargin{1},'irf.TimeTable')
  tt=varargin{1};
  if isfield(tt,'UserData')
    ud=tt.UserData;
    if isfield(ud,'directory')
      isSpecifiedDataDirectory = 1;
      handles.dirNames=ud.directory;
    end
  end
  if ~isSpecifiedDataDirectory
    disp('irf.tt_browse(): directory of time table events not specified');
    handles.dirNames=cell(1,numel(tt));
  end
  tint=tt.TimeInterval;
  handles.timeInterval = irf_time(tint,'tint>utc');
  handles.tt=tint;
  handles.userdata.figure=figure;
  handles.userdata.TTselected=irf.TimeTable;
else
  irf_log('fcal','no time table specified as input');
  return
end
handles.list_index_match_filter=1:numel(handles.dirNames);
set(handles.time_table_to_browse,'string',handles.timeInterval);
set(handles.edit1,'string','');
if numel(varargin) >= 2 % plot function specified as 2nd argument
  if ischar(varargin{2})
    set(handles.edit2,'string',varargin{2});
  end
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CAA_GUI_LIST wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tt_browse_OutputFcn(hObject, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in time_table_to_browse.
function time_table_to_browse_Callback(hObject, ~, handles)
% hObject    handle to time_table_to_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles   ~ structure with handles and user data (see GUIDATA)
i=get(handles.time_table_to_browse,'Value');
iEvent=handles.list_index_match_filter(i);
eventDir=handles.dirNames{iEvent};
tint=[handles.tt(iEvent,1) handles.tt(iEvent,2)];
cd(handles.workDirectory);
plotFunc=get(handles.edit2,'string');
[plotFuncDir,~,~] = fileparts(which(plotFunc));
if strcmp(plotFuncDir,pwd)
  addpath(pwd)
end
hcf=handles.userdata.figure; % figure handle in which plotFunc should plot
if~isempty(eventDir)
  cd(eventDir);
  run(plotFunc);
  cd(handles.workDirectory);
else
  run(plotFunc)
end
% Hints: contents = cellstr(get(hObject,'String')) returns time_table_to_browse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from time_table_to_browse


% --- Executes during object creation, after setting all properties.
function time_table_to_browse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_table_to_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filter=regexpi(get(handles.edit1,'string'),'\w*','match'); % cell array with filter values
handles.list_index_match_filter=[];
if isempty(filter)
  handles.list_index_match_filter=1:numel(handles.dirNames);
else
  for i=1:size(handles.timeInterval,1)
    ok=1;
    a=regexpi(handles.timeInterval(i,:),filter);
    for j=1:numel(a)
      ok=ok*(~isempty(a{j}));
    end
    if ok, handles.list_index_match_filter(end+1)=i; end
  end
end
if isempty(handles.list_index_match_filter)
  disp('WARNING! Filter did not match anything!');
  handles.list_index_match_filter=1:numel(handles.timeInterval);
end
set(handles.time_table_to_browse,'listboxtop',1);
set(handles.time_table_to_browse,'string',handles.timeInterval(handles.list_index_match_filter,:));
set(handles.time_table_to_browse,'value',1);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in time_table_interesting_events.
function time_table_interesting_events_Callback(hObject, eventdata, handles)
% hObject    handle to time_table_interesting_events (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iSelected=get(handles.time_table_interesting_events,'Value');
tint=handles.userdata.TTselected.TimeInterval(iSelected,:);
cd(handles.workDirectory);
plotFunc=get(handles.edit2,'string');
[plotFuncDir,~,~] = fileparts(which(plotFunc));
if strcmp(plotFuncDir,pwd)
  addpath(pwd)
end
hcf=handles.userdata.figure; % figure handle in which plotFunc should plot
run(plotFunc);

% Hints: contents = cellstr(get(hObject,'String')) returns time_table_interesting_events contents as cell array
%        contents{get(hObject,'Value')} returns selected item from time_table_interesting_events


% --- Executes during object creation, after setting all properties.
function time_table_interesting_events_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_table_interesting_events (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
    get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on button press in addevent.
function addevent_Callback(hObject, ~, handles)
% hObject    handle to addevent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hf=getfield(get(handles.userdata.figure,'userdata'),'subplot_handles');
iSelected=get(handles.addevent,'value');
switch iSelected
  case 1 % Add event
    tStart=getfield(get(handles.userdata.figure,'userdata'),'t_start_epoch');
    tint_add=tStart+get(hf(1),'xlim');
    handles.userdata.TTselected=add(handles.userdata.TTselected,tint_add);
    %handles.userdata.TTselected=unique(handles.userdata.TTselected); % if only unique events wanted
    irf_log('fcal',['Added to time table of interest: ' irf_time(tint_add,'tint>utc')]);
  case 2 % remove event
    iRemove=get(handles.time_table_interesting_events,'value');
    irf_log('fcal',['Removing interval #' num2str(iRemove) ': ' ...
      irf_time(handles.userdata.TTselected.TimeInterval(iRemove,:),'tint>utc')]);
    handles.userdata.TTselected=remove(handles.userdata.TTselected,iRemove);
    set(handles.time_table_interesting_events,'value',max(iRemove-1,1));
  case 3 % sort events
    irf_log('fcal','Sorting selected time intervals');
    handles.userdata.TTselected=sort(handles.userdata.TTselected);
end
ttsel=handles.userdata.TTselected.TimeInterval;
timeInterval_InterestingEvents=irf_time(ttsel,'tint>utc');
set(handles.time_table_interesting_events,'string',timeInterval_InterestingEvents);
irf_log('fcal',['Together ' num2str(numel(handles.userdata.TTselected)) ' events']);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in exporttable.
function exporttable_Callback(hObject, ~, handles)
% hObject    handle to exporttable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iSelected=get(handles.exporttable,'value');
switch iSelected
  case 1 % export table
    [FileName,PathName,FilterIndex] = uiputfile('*.tt','Select Time Table file (extension tt)');
    if FilterIndex
      export_ascii(handles.userdata.TTselected,[PathName FileName]);
      irf_log('dsrc',['Exporting time table to file:' PathName FileName]);
    end
  case 2 % import table
    [FileName,PathName,FilterIndex] = uigetfile('*.tt','Select Time Table file');
    if FilterIndex
      handles.userdata.TTselected=irf.TimeTable([PathName FileName]);
      irf_log('dsrc',['Imported time table from file:' PathName FileName]);
    end
    set(handles.time_table_interesting_events,'string',irf_time(handles.userdata.TTselected.TimeInterval,'tint>utc'));
    set(handles.time_table_interesting_events,'value',1);
  case 3 % assign to base
    assignin('base','TT_selected',handles.userdata.TTselected);
    disp('****************************')
    disp('Selected time intervals are')
    disp('in variable TT_selected');
    disp('****************************')
end
% Update handles structure
guidata(hObject, handles);
