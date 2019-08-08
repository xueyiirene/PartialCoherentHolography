function varargout = FastFindPeaks2D(varargin)
% FASTFINDPEAKS2D MATLAB code for FastFindPeaks2D.fig
%      FASTFINDPEAKS2D, by itself, creates a new FASTFINDPEAKS2D or raises the existing
%      singleton*.
%
%      H = FASTFINDPEAKS2D returns the handle to a new FASTFINDPEAKS2D or the handle to
%      the existing singleton*.
%
%      FASTFINDPEAKS2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FASTFINDPEAKS2D.M with the given input arguments.
%
%      FASTFINDPEAKS2D('Property','Value',...) creates a new FASTFINDPEAKS2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FastFindPeaks2D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FastFindPeaks2D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FastFindPeaks2D

% Last Modified by GUIDE v2.5 13-Jun-2019 15:36:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FastFindPeaks2D_OpeningFcn, ...
                   'gui_OutputFcn',  @FastFindPeaks2D_OutputFcn, ...
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


% --- Executes just before FastFindPeaks2D is made visible.
function FastFindPeaks2D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FastFindPeaks2D (see VARARGIN)

% Choose default command line output for FastFindPeaks2D
handles.output = hObject;
handles.Px=[];
handles.Py=[];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FastFindPeaks2D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FastFindPeaks2D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function Threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.BGthres=get(hObject,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in FindPeaks.
function FindPeaks_Callback(hObject, eventdata, handles)
% hObject    handle to FindPeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.A=double(imrotate(imread([handles.path '\' handles.ImageName '.tif']),handles.RotationAngle,'bilinear','crop')) ; % read in first image
% handles.A = eval(handles.ImageName);
axes(handles.displayImage);
cla;
imagesc(handles.A);colormap('gray');
title('choose ROI');
axis image;
% This round selection
A1=handles.A;
Caxis_auto=[];
% Autoselection spots
[rioX,rioY]=ginput(2);
rioX=round(rioX);rioY=round(rioY);
Caxis_temp=FastPeakFind(A1(rioY(1):rioY(2),rioX(1):rioX(2)),handles.BGthres*max(max(A1(rioY(1):rioY(2),rioX(1):rioX(2)))));

%pad spots
% Caxis_temp1=Caxis_temp-3;
% Caxis_temp=[Caxis_temp;Caxis_temp1];

Caxis1=[Caxis_temp(1:2:end)+rioX(1) Caxis_temp(2:2:end)+rioY(1)];
    if isempty(Caxis_auto)
       Caxis_auto=Caxis1;
    else
       Caxis_auto=[Caxis_auto;Caxis1];
    end
% Sparse spots
%     if ~isempty(Caxis_auto)
%        Caxis_auto=SparseSpots( Caxis_auto, A1, handles.radius);  
%        hold on;
%        scatter(Caxis_auto(:,1),Caxis_auto(:,2),'o');pause(0.1);
%     end
if ~isempty(Caxis_auto)
hold on;
scatter(Caxis_auto(:,1),Caxis_auto(:,2),'o');
handles.Caxis_auto=Caxis_auto;
end
guidata(hObject, handles);



function LoadImage_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LoadImage as text
%        str2double(get(hObject,'String')) returns contents of LoadImage as a double
handles.path=get(hObject,'String');
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function LoadImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveselection.
function saveselection_Callback(hObject, eventdata, handles)
% hObject    handle to saveselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Px=handles.Px+handles.Xshift;
handles.Py=handles.Py+handles.Yshift;
h=round([handles.Px/handles.Xratio handles.Py/handles.Yratio]);
uisave('h',[handles.path '\FindPeaks.mat']);
guidata(hObject, handles);


% --- Executes on button press in keepselection.
function keepselection_Callback(hObject, eventdata, handles)
% hObject    handle to keepselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.Px)
handles.Py=cat(1,handles.Py,handles.Caxis_auto(:,1));
handles.Px=cat(1,handles.Px,handles.Caxis_auto(:,2));
else
    handles.Py=handles.Caxis_auto(:,1);
    handles.Px=handles.Caxis_auto(:,2);
end
guidata(hObject,handles);



function Xratio_Callback(hObject, eventdata, handles)
% hObject    handle to Xratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xratio as text
%        str2double(get(hObject,'String')) returns contents of Xratio as a double
handles.Xratio=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Xratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ImageName_Callback(hObject, eventdata, handles)
% hObject    handle to ImageName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ImageName=get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ImageName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RotationAngle_Callback(hObject, eventdata, handles)
% hObject    handle to RotationAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.RotationAngle=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function RotationAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RotationAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Yratio_Callback(hObject, eventdata, handles)
% hObject    handle to Yratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Yratio=str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Yratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Yratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Xshift_Callback(hObject, eventdata, handles)
% hObject    handle to Xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Xshift=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Xshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Yshift_Callback(hObject, eventdata, handles)
% hObject    handle to Yshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Yshift=str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Yshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Yshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
