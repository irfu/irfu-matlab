function varargout = caa_gui_list(varargin)
% CAA_GUI_LIST list/filter interactively cell array values
% CAA_GUI_LIST(list,values)
%      CAA_GUI_LIST, by itself, creates a new CAA_GUI_LIST or raises the existing
%      singleton*.
%
%      H = CAA_GUI_LIST returns the handle to a new CAA_GUI_LIST or the handle to
%      the existing singleton*.
%
%      CAA_GUI_LIST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAA_GUI_LIST.M with the given input arguments.
%
%      CAA_GUI_LIST('Property','Value',...) creates a new CAA_GUI_LIST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CAA_GUI_LIST_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CAA_GUI_LIST_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CAA_GUI_LIST

% Last Modified by GUIDE v2.5 28-Jun-2012 14:02:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @caa_gui_list_OpeningFcn, ...
                   'gui_OutputFcn',  @caa_gui_list_OutputFcn, ...
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


% --- Executes just before caa_gui_list is made visible.
function caa_gui_list_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to caa_gui_list (see VARARGIN)

% Choose default command line output for CAA_GUI_LIST
handles.output = hObject;
handles.list=varargin{1};
handles.list_index_match_filter=1:numel(handles.list);
handles.values=varargin{2};
set(handles.listbox1,'string',handles.list);
set(handles.listbox2,'string',handles.values(:,1));
set(handles.edit1,'string','');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CAA_GUI_LIST wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = caa_gui_list_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=get(handles.listbox1,'Value');
set(handles.listbox2,'value',i);
set(handles.edit2,'string',edit2_output(handles));

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filter=regexpi(get(handles.edit1,'string'),'\w*','match'); % cell array with filter values
handles.list_index_match_filter=[];
if isempty(filter)
	handles.list_index_match_filter=1:numel(handles.list);
else
	for i=1:numel(handles.list)
		ok=1;
		a=regexpi(handles.list{i},filter);
		for j=1:numel(a)
			ok=ok*(~isempty(a{j}));
		end
		if ok, handles.list_index_match_filter(end+1)=i; end
	end
end
if isempty(handles.list_index_match_filter)
	disp('WARNING! Filter did not match anything!');
	handles.list_index_match_filter=1:numel(handles.list);
end
set(handles.listbox1,'listboxtop',1);
set(handles.listbox1,'string',handles.list(handles.list_index_match_filter));
set(handles.listbox1,'value',1);
set(handles.listbox2,'listboxtop',1);
set(handles.listbox2,'string',handles.values(handles.list_index_match_filter));
set(handles.listbox2,'value',1);
set(handles.edit2,'string',edit2_output(handles));
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

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


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
i=get(handles.listbox2,'Value');
set(handles.listbox1,'value',i);
set(handles.edit2,'string',edit2_output(handles));

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
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

function txt=edit2_output(handles)
% return what to put in edit2 text field
i=get(handles.listbox1,'Value');
ii=handles.list_index_match_filter(i);
txt={};
txt{1}=handles.list{ii};
for j=1:size(handles.values,2)
	txt{end+1}='';
	txt{end+1}=handles.values{ii,j};
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
