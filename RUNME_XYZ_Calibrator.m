function varargout = RUNME_XYZ_Calibrator(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @RUNME_XYZ_Calibrator_OpeningFcn, ...
    'gui_OutputFcn',  @RUNME_XYZ_Calibrator_OutputFcn, ...
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
end

function RUNME_XYZ_Calibrator_OpeningFcn(hObject, eventdata, handles, varargin)
clc
disp('Welcome to Partially Coherenet Holography Calibration GUI')
disp('Nidaq Starting')
handles.Setup = Function_Load_Parameters();
disp('Started')
global DAQstate;DAQstate=[0 0 0 0 0 0];
global status; status = 1; % 1 if ok to continue
global ToSave;
global Stage;
global Cam; Cam.State=0;
global StageStatus; StageStatus=0;
global stageposition; stageposition= [0 0 0];
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider1,'Value',0); %Reset slider
set(handles.slider2,'Value',0); % Reset slider
set(handles.edit1,'string',0); %Initial stim laser Voltage
set(handles.edit2,'string',0); %Initial Imaging laser voltageset(handles.edit6,'string',int2str(10)); %Number of pulses
set(handles.edit3,'string','OFF');%Stage X
set(handles.edit4,'string','OFF');%Stage Y
set(handles.edit5,'string','OFF');%Stage Z
set(handles.edit6,'string','Camera OFF');%Camera status field
set(handles.edit7,'string',0); %DMD Target X location in pixels
set(handles.edit8,'string',0); %DMD Target Y location in pixels
set(handles.edit9,'string',0); %DMD Target Z location in pixels Not defined in the beginign.
set(handles.edit10,'string',12);%DMD Target radius location in pixels
set(handles.edit11,'string',100);%Range in microns for substage calibration
set(handles.edit12,'string',8);%Number of Points for substage calibration
set(handles.edit13,'string',1.7);%Number of camera pixels per microns
set(handles.edit14,'string',1.6);%Voltage for DMD cam calibraition
set(handles.edit15,'string',150);%DMD range pixels for DMD cam calibraition
set(handles.edit16,'string',3);%Number of pixels for DMD cam calibraition
set(handles.edit17,'string',30);%DeltaZ for spatial DMD calibration
set(handles.edit18,'string',3);%Max Galvo voltage (volts)
set(handles.edit19,'string',3);%Number of voltage steps at galvos
Stim.BlankFrame = ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY);
function_directfeed_DMD( handles.Setup,Stim.BlankFrame );
handles.output = hObject;
guidata(hObject, handles);
end


% --- Executes DMD GalvoCalibration
function pushbutton18_Callback(hObject, eventdata, handles)
global Cam;
global status;status = 1;
global StageStatus;
global DAQstate;
global Stage;
global stageposition;
axes(handles.axes1); hold off; cla;axes(handles.axes2); hold off; cla;
DAQstate(1) = str2num(get(handles.edit14,'string'));
axes(handles.axes1); hold off; cla;axes(handles.axes2); hold off; cla;
    UX = linspace(0,(handles.Setup.DMD.LX-1),handles.Setup.DMD.LX);
    UY = linspace(0,(handles.Setup.DMD.LY-1),handles.Setup.DMD.LY);
    UX = UX-mean(UX); UY = UY-mean(UY);
    [XX,YY] = ndgrid(UX,UY);
    Frame = double((XX).^2 +(YY).^2 < (str2double(get(handles.edit10,'string'))).^2);
    function_directfeed_DMD( handles.Setup,Frame ); outputSingleScan(handles.Setup.Daq,DAQstate);

Uteta = linspace(-str2double(get(handles.edit18,'string')),str2double(get(handles.edit18,'string')),floor(str2double(get(handles.edit19,'string'))));
[UTX,UTY] = ndgrid(Uteta,Uteta);
Galvo.VoltageXYPts = [UTX(:) UTY(:)];
[LN,~] = size(Galvo.VoltageXYPts);
Galvo.CAMPts = Galvo.VoltageXYPts-Galvo.VoltageXYPts;
axes(handles.axes2); scatter(Galvo.VoltageXYPts(:,1), Galvo.VoltageXYPts(:,2)); xlabel('Galvo x'); ylabel('Galvo Y');hold on ; pause(0.5)


if Cam.State == 1;
    if StageStatus ==1;
    stageposition = [0 0 0];
            stageposition(3) = str2num(get(handles.edit17,'string'));
            set(handles.edit3,'string',num2str(stageposition(1)));%Stage X
            set(handles.edit4,'string',num2str(stageposition(2)));%Stage Y
            set(handles.edit5,'string',num2str(stageposition(3)));%Stage Z
            function_Goto_stage( Stage,stageposition ); 
            
            
            
        
    for i = 1:LN
DAQstate(3) = Galvo.VoltageXYPts(i,1);
DAQstate(4) = Galvo.VoltageXYPts(i,2);
outputSingleScan(handles.Setup.Daq,DAQstate);
[ Data ] = function_GetFrameCam( Cam.cam,Cam.Parameters,1 );
            axes(handles.axes2); scatter( Galvo.VoltageXYPts(i,1),  Galvo.VoltageXYPts(i,2),'filled','red'); hold on
            axes(handles.axes1); imagesc(Data(:,:,1));  colormap gray; caxis([0 255]); hold on;
            pause(0.05)
            %imagesc(Frame);pause(1)
            if status ==0;  set(handles.edit6,'string','Camera ON'); disp('Interrupted'); return; end;
            [ x,y ] = function_find_blob( squeeze(Data(:,:,1)),10);
            scatter(y,x,'red');title(int2str(i)); %%%%%% VERIFY THIS
            pause(0.1);
            Galvo.CAMPts(i,1) = x;
            Galvo.CAMPts(i,2) = y;
            hold off;
    
    end
    
    stageposition = [0 0 0];        
            set(handles.edit3,'string',num2str(stageposition(1)));%Stage X
            set(handles.edit4,'string',num2str(stageposition(2)));%Stage Y
            set(handles.edit5,'string',num2str(stageposition(3)));%Stage Z
            function_Goto_stage( Stage,stageposition ); 
            
    DAQstate = [0 0 0 0 0 0]; outputSingleScan(handles.Setup.Daq,DAQstate);
%%% Proper ocnversion Syntax... 
Guessgalvo = function_CofC( Galvo.CAMPts, Galvo.VoltageXYPts, Galvo.CAMPts );
GuessCam = function_CofC(  Galvo.VoltageXYPts, Galvo.CAMPts, Galvo.VoltageXYPts);

axes(handles.axes1); hold off; cla; 
scatter(Galvo.VoltageXYPts(:,1), Galvo.VoltageXYPts(:,2),'red'); hold on;
scatter(Guessgalvo(:,1), Guessgalvo(:,2),'blue','filled'); hold off;
axes(handles.axes2); hold off; cla;
scatter(Galvo.CAMPts(:,1), Galvo.CAMPts(:,2),'red'); hold on
scatter(GuessCam(:,1), GuessCam(:,2),'blue','filled'); hold off;        
            
            
            
            
            
    else
        set(handles.edit6,'string','Turn Stage On first !');%Camera status field
    end      
else
    set(handles.edit6,'string','Turn Camera On first !');%Camera status field
end


end

% --- Executes DMD CAM Calibration
function pushbutton17_Callback(hObject, eventdata, handles)
global Cam;
global status;status = 1; 
global DAQstate;
axes(handles.axes1); hold off; cla;axes(handles.axes2); hold off; cla;
DAQstate(1) = str2num(get(handles.edit14,'string'));
axes(handles.axes1); hold off; cla;axes(handles.axes2); hold off; cla;
    UX = linspace(0,(handles.Setup.DMD.LX-1),handles.Setup.DMD.LX);
    UY = linspace(0,(handles.Setup.DMD.LY-1),handles.Setup.DMD.LY);
    UX = UX-mean(UX); UY = UY-mean(UY);
    [XX,YY] = ndgrid(UX,UY);
   

if Cam.State == 1;
% genreate maze of DMD points, and match them to Cam points
r = linspace(-str2num(get(handles.edit15,'string')),str2num(get(handles.edit15,'string')),floor(str2num(get(handles.edit16,'string'))));
[UX,UY] = ndgrid(r,r);
Points.DMD = [UX(:) UY(:)];
[LN,~] = size(Points.DMD);
Points.CAM = Points.DMD-Points.DMD;
axes(handles.axes2); scatter(Points.DMD(:,1), Points.DMD(:,2)); hold on
for j = 1:LN
  
    Frame = double((XX-Points.DMD(j,1)).^2 +(YY-Points.DMD(j,2)).^2 < (str2double(get(handles.edit10,'string'))).^2);
    function_directfeed_DMD( handles.Setup,Frame ); outputSingleScan(handles.Setup.Daq,DAQstate);
            [ Data ] = function_GetFrameCam( Cam.cam,Cam.Parameters,1 );
            axes(handles.axes2); scatter(Points.DMD(j,1), Points.DMD(j,2),'filled','red'); hold on
            axes(handles.axes1); imagesc(Data(:,:,1));  colormap gray; caxis([0 255]); hold on;
            pause(0.05)
            %imagesc(Frame);pause(1)
            if status ==0;  set(handles.edit6,'string','Camera ON'); disp('Interrupted'); return; end;
            [ x,y ] = function_find_blob( squeeze(Data(:,:,1)),10);
            scatter(y,x,'red');title(int2str(j)); %%%%%% VERIFY THIS
            Points.CAM(j,1) = x;
            Points.CAM(j,2) = y;
            hold off;
    pause(0.1);
    
end


DAQstate = [0 0 0 0 0 0]; outputSingleScan(handles.Setup.Daq,DAQstate);
%%% Proper ocnversion Syntax... 
GuessCam = function_CofC( Points.DMD, Points.CAM, Points.DMD );
GuessDMD = function_CofC( Points.CAM, Points.DMD, Points.CAM );

axes(handles.axes1); hold off; cla; 
scatter(Points.DMD(:,1), Points.DMD(:,2),'red'); hold on;
scatter(GuessDMD(:,1), GuessDMD(:,2),'blue','filled'); hold off;
axes(handles.axes2); hold off; cla;
scatter(Points.CAM(:,1), Points.CAM(:,2),'red'); hold on
scatter(GuessCam(:,1), GuessCam(:,2),'blue','filled'); hold off;

try 
    load('GUI_DMD_Camera_Calibraiton.mat');
catch;    
end;
SaveCalib.DMDCAMPoints = Points;    
save('GUI_DMD_Camera_Calibraiton.mat','SaveCalib');       

   
else
    set(handles.edit6,'string','Turn Camera On first !');%Camera status field
end

end



% --- Executes Substage_Calibration
function pushbutton16_Callback(hObject, eventdata, handles)
global status; status = 1;
global Cam;
global Stage;
global DAQstate;
global StageStatus;
global stageposition;
 UX = linspace(0,(handles.Setup.DMD.LX-1),handles.Setup.DMD.LX);
    UY = linspace(0,(handles.Setup.DMD.LY-1),handles.Setup.DMD.LY);
    UX = UX-mean(UX); UY = UY-mean(UY);
    [XX,YY] = ndgrid(UX,UY);
DAQstate(1) = str2num(get(handles.edit14,'string'));
Frame = double((XX).^2 +(YY).^2 < (str2double(get(handles.edit10,'string'))).^2);
    function_directfeed_DMD( handles.Setup,Frame ); outputSingleScan(handles.Setup.Daq,DAQstate);
SubstagePoints = function_cross(str2double(get(handles.edit3,'string')),str2double(get(handles.edit4,'string')),str2double(get(handles.edit11,'string')),str2double(get(handles.edit11,'string')),str2double(get(handles.edit12,'string')));
[~,LP] = size(SubstagePoints);
CamPoints = SubstagePoints-SubstagePoints;
axes(handles.axes1); hold off; cla;axes(handles.axes2); hold off; cla;
scatter(SubstagePoints(1,:),SubstagePoints(2,:),'filled','red'); xlabel('Substage dX um'); ylabel('Substage dY um'); pause(1);
if Cam.State == 1;
    if StageStatus ==1;
        set(handles.edit6,'string','Substage calibration');%Camera status field
    %Set to zero
        function_Zero_stage( Stage );
    stageposition = [0 0 0];
    set(handles.edit3,'string',num2str(0));%Stage X
    set(handles.edit4,'string',num2str(0));%Stage Y
    set(handles.edit5,'string',num2str(0));%Stage Z
    
        for j = 1:LP
            stageposition = [0 0 0];
            stageposition(1) = SubstagePoints(1,j);
            stageposition(2) = SubstagePoints(2,j);
            set(handles.edit3,'string',num2str(stageposition(1)));%Stage X
            set(handles.edit4,'string',num2str(stageposition(2)));%Stage Y
            set(handles.edit5,'string',num2str(stageposition(3)));%Stage Z
            function_Goto_stage( Stage,stageposition ); 
            pause(1);
            [ Data ] = function_GetFrameCam( Cam.cam,Cam.Parameters,1 );
            axes(handles.axes1);
            imagesc(Data(:,:,1));  colormap gray; caxis([0 255]); hold on
            if status ==0;  set(handles.edit6,'string','Camera ON'); return; end;
            [ x,y ] = function_find_blob( squeeze(Data(:,:,1)),10);
            scatter(y,x,'red'); %%%%%% VERIFY THIS
            CamPoints(1,j) = x;
            CamPoints(2,j) = y;
            hold off;
            axes(handles.axes2);
            scatter(y,x,'red','filled'); hold on;
        end
        stageposition = [0 0 0];
        set(handles.edit3,'string',num2str(stageposition(1)));%Stage X
        set(handles.edit4,'string',num2str(stageposition(2)));%Stage Y
        set(handles.edit5,'string',num2str(stageposition(3)));%Stage Z
        function_Goto_stage( Stage,stageposition ); pause(0.2);
 
        %%%%%% End processign the data here... 
        CamPoints(1,:) = CamPoints(1,:)-CamPoints(1,1);
        CamPoints(2,:) = CamPoints(2,:)-CamPoints(2,1);
        CamRadius = sqrt(CamPoints(1,:).^2+CamPoints(2,:).^2);
        SubstagePoints(1,:)=SubstagePoints(1,:)-SubstagePoints(1,1);
        SubstagePoints(2,:) =SubstagePoints(2,:)-SubstagePoints(1,1);
        StageRadius = sqrt( SubstagePoints(1,:).^2+ SubstagePoints(2,:).^2);
        keep = StageRadius>20;
        CamRadius = CamRadius(keep);
        StageRadius = StageRadius(keep);
        FirstCalibration.pixelpermicrons = mean(CamRadius./StageRadius);
        set(handles.edit13,'string',FirstCalibration.pixelpermicrons);
        axes(handles.axes1); hold off; cla; scatter(CamRadius,StageRadius); xlabel('Camera D-Pixels'),ylabel('Stage displacement')
        
        try 
    load('GUI_DMD_Camera_Calibraiton.mat');
        catch;
        end;
SaveCalib.CamPixelPermicrons = FirstCalibration.pixelpermicrons;    
save('GUI_DMD_Camera_Calibraiton.mat','SaveCalib');       

        
        
    else
        set(handles.edit6,'string','Turn Stage On first !');%Camera status field
    end
else
    set(handles.edit6,'string','Turn Camera On first !');%Camera status field
end


end



% --- Executes DMD OPEN
function pushbutton13_Callback(hObject, eventdata, handles)
Frame = ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY);
function_directfeed_DMD( handles.Setup,Frame );
global DAQstate;
outputSingleScan(handles.Setup.Daq,DAQstate);
axes(handles.axes1); hold off; cla; imagesc(Frame); caxis([0 1]);
end
% --- Executes DMD Close
function pushbutton14_Callback(hObject, eventdata, handles)
Frame = zeros(handles.Setup.DMD.LX,handles.Setup.DMD.LY);
function_directfeed_DMD( handles.Setup,Frame );
global DAQstate;
outputSingleScan(handles.Setup.Daq,DAQstate);
axes(handles.axes1); hold off; cla; imagesc(Frame); caxis([0 1]);
end
% --- Executes DMD Target
function pushbutton15_Callback(hObject, eventdata, handles)
axes(handles.axes1); hold off; cla;
UX = linspace(0,(handles.Setup.DMD.LX-1),handles.Setup.DMD.LX);
UY = linspace(0,(handles.Setup.DMD.LY-1),handles.Setup.DMD.LY);
UX = UX-mean(UX); UY = UY-mean(UY);
[XX,YY] = ndgrid(UX,UY);
Frame = double((XX-str2double(get(handles.edit7,'string'))).^2 +(YY-str2double(get(handles.edit8,'string'))).^2 < (str2double(get(handles.edit10,'string'))).^2);
imagesc(Frame); caxis([0 1]);
function_directfeed_DMD( handles.Setup,Frame );
global DAQstate;
outputSingleScan(handles.Setup.Daq,DAQstate);
end


% --- Executes Start Cam
function pushbutton10_Callback(hObject, eventdata, handles)
global Cam;
global status; status = 1;
try
   if Cam.State == 0
    [ Cam.cam,Cam.Parameters ] = function_StartCam( );
    Cam.State = 1;
    set(handles.edit6,'string','Camera ON');%Camera status field
   elseif Cam.State == 1
     set(handles.edit6,'string','Camera Already ON');%Camera status field  
   end
catch
    set(handles.edit6,'string','Camera Startup Error');%Camera status field
    Cam.State = -1;
end
end

% --- Executes Stop Cam
function pushbutton11_Callback(hObject, eventdata, handles)
global Cam;
try
    if Cam.State ==0; set(handles.edit6,'string','Camera Already OFF');%Camera status field
    elseif Cam.State ==1;
        Cam.cam.Exit; Cam.State =0;
        set(handles.edit6,'string','Camera OFF');
    else
    end
catch
    set(handles.edit6,'string','Camera Startup Error');%Camera status field
    Cam.State = -1;
end
end

% --- Executes Live Preview
function pushbutton12_Callback(hObject, eventdata, handles)
global status; status = 1;
global Cam;
axes(handles.axes1); hold off; cla; axes(handles.axes2); hold off; cla;
if Cam.State == 1;
    set(handles.edit6,'string','Camera Preview');%Camera status field
    Cam.State=2;
    for j = 1:1000000
        [ Data ] = function_GetFrameCam( Cam.cam,Cam.Parameters,1 );
        axes(handles.axes1);
        imagesc(Data(:,:,1));  colormap gray; caxis([0 255]); title(['Frame' int2str(j)]);
        axes(handles.axes2);
        imagesc(Data((handles.Setup.Cam.LX/2-150):(handles.Setup.Cam.LX/2+150),(handles.Setup.Cam.LY/2-150):(handles.Setup.Cam.LY/2+150),1));  colormap gray; caxis([0 255]); title(['Zoom']);
        pause(0.1);
        if status ==0;  set(handles.edit6,'string','Camera ON'); Cam.State=1; return; end;
    end
else
    set(handles.edit6,'string','Turn camera on first !');%Camera status field
end
end

function varargout = RUNME_XYZ_Calibrator_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
end


% --- Executes Slider STIM laser
function slider1_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
set(handles.edit1,'string',round(10*h)/10);
global DAQstate;
DAQstate(1) = h;
outputSingleScan(handles.Setup.Daq,DAQstate);
end

% --- Executes Slider imaging laser
function slider2_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
set(handles.edit2,'string',round(10*h)/10);
global DAQstate;
DAQstate(2) = h;
outputSingleScan(handles.Setup.Daq,DAQstate);
end


% --- Executes Lasers off.
function pushbutton1_Callback(hObject, eventdata, handles)
global DAQstate;
DAQstate = [0 0 0 0 0 0];
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider1,'Value',0)
set(handles.slider2,'Value',0)
set(handles.edit1,'string',0);
set(handles.edit2,'string',0);
end

% --- Executes Stage ON
function pushbutton2_Callback(hObject, eventdata, handles)
global Stage;
global StageStatus;
try
    Stage  = function_Start_stage( handles.Setup.MechStageComport );
    function_Zero_stage( Stage );
    StageStatus = 1;
    set(handles.edit3,'string',num2str(0));%Stage X
    set(handles.edit4,'string',num2str(0));%Stage Y
    set(handles.edit5,'string',num2str(0));%Stage Z
catch
    set(handles.edit3,'string','ERR');%Stage X
    set(handles.edit4,'string','ERR');%Stage Y
    set(handles.edit5,'string','ERR');%Stage Z
end
end

% --- Executes Stage OFF.
function pushbutton3_Callback(hObject, eventdata, handles)
global Stage;
global StageStatus;
try
    StageStatus=0;
    function_Stop_stage( Stage );
    set(handles.edit3,'string','OFF');%Stage X
    set(handles.edit4,'string','OFF');%Stage Y
    set(handles.edit5,'string','OFF');%Stage Z
catch
    set(handles.edit3,'string','ERR');%Stage X
    set(handles.edit4,'string','ERR');%Stage Y
    set(handles.edit5,'string','ERR');%Stage Z
end
end

% --- Executes Set stage zero.
function pushbutton4_Callback(hObject, eventdata, handles)
global Stage;
global stageposition;
global StageStatus;
try
    function_Zero_stage( Stage );
    stageposition = [0 0 0];
    set(handles.edit3,'string',num2str(0));%Stage X
    set(handles.edit4,'string',num2str(0));%Stage Y
    set(handles.edit5,'string',num2str(0));%Stage Z
catch
    set(handles.edit3,'string','ERR');%Stage X
    set(handles.edit4,'string','ERR');%Stage Y
    set(handles.edit5,'string','ERR');%Stage Z
end
end

% --- Executes Read stage position.
function pushbutton5_Callback(hObject, eventdata, handles)
try
    global Stage;
    global stageposition;
    [ positionmicron ] = function_Get_stage(Stage );
    stageposition = positionmicron;
    set(handles.edit3,'string',num2str(positionmicron(1)));%Stage X
    set(handles.edit4,'string',num2str(positionmicron(2)));%Stage Y
    set(handles.edit5,'string',num2str(positionmicron(3)));%Stage Z
catch
    set(handles.edit3,'string','ERR');%Stage X
    set(handles.edit4,'string','ERR');%Stage Y
    set(handles.edit5,'string','ERR');%Stage Z
end
end

% --- Executes Go to zero
function pushbutton6_Callback(hObject, eventdata, handles)
try
    global Stage;
    position = [0 0 0]; function_Goto_stage( Stage,position ); pause(0.2);
catch
    set(handles.edit3,'string','ERR');%Stage X
    set(handles.edit4,'string','ERR');%Stage Y
    set(handles.edit5,'string','ERR');%Stage Z
end
end


% --- Executes Go to position
function pushbutton7_Callback(hObject, eventdata, handles)
try
    global Stage;
    global stageposition;
    stageposition(1) = str2num(get(handles.edit3,'string'));
    stageposition(2) = str2num(get(handles.edit4,'string'));
    stageposition(3) = str2num(get(handles.edit5,'string'));
    function_Goto_stage( Stage,stageposition ); pause(0.2);
catch
    set(handles.edit3,'string','ERR');%Stage X
    set(handles.edit4,'string','ERR');%Stage Y
    set(handles.edit5,'string','ERR');%Stage Z
end
end

% --- Executes Quit
function pushbutton8_Callback(hObject, eventdata, handles)
%Lasers off, stage Off, camera Off
global Stage;
global Cam;
global DAQstate;
DAQstate = [0 0 0 0 0 0];
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider1,'Value',0)
set(handles.slider2,'Value',0)
set(handles.edit1,'string',0);
set(handles.edit2,'string',0);
try
    function_Goto_stage( Stage,[ 0 0 0] ); pause(0.5);
    function_Stop_stage( Stage );
catch
end

try
    Cam.cam.Exit;
catch
end
close;
end

% --- Executes Interrupt
function pushbutton9_Callback(hObject, eventdata, handles)
global status;
status = 0; % 1 if ok to continue
end

%USELESS stuff down below
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end
function edit1_Callback(hObject, eventdata, handles)
end
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit2_Callback(hObject, eventdata, handles)
end
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit3_Callback(hObject, eventdata, handles)
end
function edit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit4_Callback(hObject, eventdata, handles)
end
function edit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit5_Callback(hObject, eventdata, handles)
end
function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit6_Callback(hObject, eventdata, handles)
end
function edit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit7_Callback(hObject, eventdata, handles)
end
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit8_Callback(hObject, eventdata, handles)
end
function edit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit9_Callback(hObject, eventdata, handles)
end
function edit9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit10_Callback(hObject, eventdata, handles)
end
function edit10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit11_Callback(hObject, eventdata, handles)
end
function edit11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit12_Callback(hObject, eventdata, handles)
end
function edit12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit13_Callback(hObject, eventdata, handles)
end
function edit13_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit14_Callback(hObject, eventdata, handles)
end
function edit14_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit15_Callback(hObject, eventdata, handles)
end
function edit15_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit16_Callback(hObject, eventdata, handles)
end
function edit16_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit17_Callback(hObject, eventdata, handles)
end
function edit17_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit18_Callback(hObject, eventdata, handles)
end
function edit18_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit19_Callback(hObject, eventdata, handles)
end
function edit19_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end