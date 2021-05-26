function varargout = Controller_v0(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Controller_v0_OpeningFcn, ...
    'gui_OutputFcn',  @Controller_v0_OutputFcn, ...
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

% --- Executes just before Controller_v0 is made visible.
function Controller_v0_OpeningFcn(hObject, eventdata, handles, varargin)
clc;
disp('Welcome to 3D-MAP GUI v0');disp('Nidaq Starting');
handles.Setup = Function_Load_Parameters('EphyS ON');
disp('Started')
global DAQstate;
global status; status = 1; % 1 if ok to continue
global SaveID; SaveID=0;
global ToSave;
global Stage;
global Cloud;
global StageStatus; StageStatus=0;
DAQstate=handles.Setup.DAQstateZero;
handles.stageposition = [0 0 0];
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0); %Reset slider
set(handles.slider2,'Value',0); % Reset slider
set(handles.edit3,'string',0); %Initial stim laser Voltage
set(handles.edit4,'string',0); %Initial Imaging laser voltageset(handles.edit6,'string',int2str(10)); %Number of pulses
set(handles.edit6,'string',int2str(10)); % Number of pulses
set(handles.edit7,'string',int2str(20)); % Frequency in Hertz
set(handles.edit8,'string',int2str(4)); % Pulse duration in ms
set(handles.edit9,'string',int2str(10)); % Number of repetitions
set(handles.edit10,'string',int2str(0)); % Delay between repeitions in seconds
set(handles.edit11,'string',num2str(1.6)); %Stim laser voltage,center (default)
set(handles.edit29,'string',num2str(2)); %Pre stimulation delay in seconds
set(handles.edit30,'string',num2str(0.8)); %Pulsed lase Duty cycle between 0 and 1
set(handles.edit12,'string',num2str(0.6)); %Voltage sweep start voltage
set(handles.edit13,'string',num2str(5));%Voltage sweep End voltage
set(handles.edit14,'string',num2str(10));%Voltage sweep number of steps
set(handles.edit16,'string',num2str(0)); %Targets positions pixels X
set(handles.edit17,'string',num2str(0)); %Targets positions pixels Y
set(handles.edit18,'string',num2str(0)); %Targets positions pixels Z
set(handles.edit19,'string',num2str(10)); %Target radius
set(handles.edit20,'string',num2str(0)); % PPSFRange of measurements on X axis in pixels,start
set(handles.edit21,'string',num2str(0));% PPSF Range of measurements on Y axis in pixels,start
set(handles.edit22,'string',num2str(0));% PPSFRange of measurements on Z axis in pixels,start
set(handles.edit47,'string',num2str(0)); % PPSFRange of measurements on X axis in pixels,end
set(handles.edit48,'string',num2str(0));% PPSF Range of measurements on Y axis in pixels,end
set(handles.edit49,'string',num2str(0));% PPSFRange of measurements on Z axis in pixels,end
set(handles.edit23,'string',num2str(0));% PPSF number of points on X axis
set(handles.edit24,'string',num2str(0));% PPSF number of points on Y axis
set(handles.edit25,'string',num2str(1));% PPSF number of points on Z axis
set(handles.edit26,'string','OFF');%Stage X
set(handles.edit27,'string','OFF');%Stage Y
set(handles.edit28,'string','OFF');%Stage Z
set(handles.edit39,'string',num2str(0.05));%Timesweep Min stim duration [ms]
set(handles.edit40,'string',num2str(3));%Timesweep Max stim duration [ms]
set(handles.edit41,'string',num2str(15));%Timesweep number of steps
set(handles.edit55,'string',num2str(2.5));% Sampling grid
set(handles.edit56,'string',num2str(5));%number of simultaneous foci in 3D
set(handles.edit57,'string',num2str(1));%number of simultaneous foci in 2D
set(handles.edit58,'string',num2str(100));%camera exposure time for preview
set(handles.edit60,'string',num2str(1));% number of images of 1 trigger
set(handles.edit61,'string',num2str(0.57));% select ROIs thresholding for adaptive imbinarize
set(handles.edit63,'string',num2str(2)); % Ca2+ laser voltage. Yellow laser voltage used only in Ca2+ imaging
set(handles.edit64,'string',num2str(3)); % Big regions to be permuted in ROI illumination function
set(handles.edit65,'string',num2str(1)); % zoom in times of preview
set(handles.edit66,'string',num2str(3.3)); % voltage for current injection
set(handles.edit67,'string',num2str(4)); % number of interleave patterns
Cloud.divider = handles.Setup.PointCloud.divider;% number of legs
Cloud.AngleMagnitude = handles.Setup.PointCloud.AngleMagnitude;% galvo voltage
Cloud.AngleMagnitudeBackup = Cloud.AngleMagnitude;
Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;
Cloud.CycleLength = handles.Setup.PointCloud.CycleLength;
set(handles.edit42,'string',int2str(handles.Setup.PointCloud.divider));% number of legs
set(handles.edit43,'string',num2str(handles.Setup.PointCloud.AngleMagnitude));% galvo voltage
set(handles.edit44,'string',num2str(handles.Setup.PointCloud.CycleLength*1000));%cycle length ms
set(handles.edit31,'string',['Dummy_File_' int2str(floor(rand(1,1)*10000))]);%Stage Z
Stim.BlankFrame = ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY);
set(handles.checkbox1,'Value', 1);
set(handles.checkbox2,'Value', 0);
set(handles.checkbox8,'Value', 0);
set(handles.checkbox11,'Value', 0);
set(handles.checkbox13,'Value', 0);
set(handles.checkbox16,'Value', 0);
set(handles.checkbox17,'Value', 0);
set(handles.checkbox18,'Value', 0);
function_directfeed_DMD( handles.Setup,Stim.BlankFrame);
handles.output = hObject;
guidata(hObject, handles);
end


% --- Executes Temporal Optimization
function pushbutton21_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);set(handles.edit4,'string',0);
Stim.Npeaks = str2num(get(handles.edit6,'string')); % Number of blue light pulses in test
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
if get(handles.checkbox3,'Value') == 1
    Stim.Timeramp = exp(linspace(log(str2num(get(handles.edit39,'string'))/1000),log(str2num(get(handles.edit40,'string'))/1000),floor(str2num(get(handles.edit41,'string')))));
else
    Stim.Timeramp = linspace(str2num(get(handles.edit39,'string'))/1000,str2num(get(handles.edit40,'string'))/1000,floor(str2num(get(handles.edit41,'string')))); %Voltage ramp to test how much light is needed to stim...
end
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
% while Stim.Output(find(Stim.UUT>Stim.DelayBorder,1,'first'),1)~=0
%     Stim.DelayBorder=Stim.DelayBorder+0.0001;
% end 
Stim.Baseline = Stim.Output(:,1);

axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
Stim.Output(end,:)=0;
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
end
Result.DMDFrames =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
function_feed_DMD(handles.Setup, uint8(gather(Result.DMDFrames))*255);

for i = 1:floor(str2num(get(handles.edit9,'string'))) % Do many repetitions
    for j = 1:numel(Stim.Timeramp)
        if status == 0; disp('Procedure interrupted');  break; end;
        Stim.Output(1:end-1,1)=Stim.Voltage;
        Stim.nonzerovalues = find(Stim.Output(:,1));
        Stim.Output(1:end-1,6)=1;
        Stim.Output(:,2:5)=0;
        Stim.Output(:,7)=0;
        for n = 1:Stim.Npeaks
            Stim.Output(:,7) = Stim.Output(:,7)+double(Stim.UUT>Stim.DelayBorder+(n-1)/Stim.FreqHZ)'.*double(Stim.UUT<Stim.DelayBorder+Stim.Timeramp(j)+(n-1)/Stim.FreqHZ)';
        end
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        Stim.baseline=mean(Data(1:100));
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
        Cell.Timescore{j,i} = score;
        DataSave{j,i} = Data;
        axes(handles.axes2); 
        scatter(Stim.Timeramp(j),Cell.Timescore{j,i},'filled','red'); hold on;xlabel('Pulse duration'); ylabel('Ephys score [AU]');title(['Calibration pass ' int2str(i) 'of' int2str(floor(str2num(get(handles.edit9,'string'))))])
        axes(handles.axes3);
        plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage'); title([ 'Score = ' num2str(floor(Cell.Timescore{j,i}*100)/100)]);
        axes(handles.axes4);
        plot(Stim.UUT,Stim.Output(:,7));xlabel('Time [s]'); ylabel('Stim laser [V]');title(['Sequence ' int2str(i) ' of ' int2str(floor(str2num(get(handles.edit9,'string'))))]);

    end
end
handles.Setup=function_StopProj_DMD(handles.Setup);
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
try
    ToSave.type = 'Temporal Optimization';
    ToSave.Stim = Stim;
    ToSave.status = status;
    ToSave.Data = DataSave;
    ToSave.VoltageScores = Cell.Timescore;
    ToSave.VoltageRamp = Stim.Voltageramp;
catch;
end;
    if get(handles.checkbox1,'Value') == 1
        filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
        save(filename,'ToSave');
        disp('Data Saved')
    end  

end

% --- Executes ForceGalvoStop
function checkbox4_Callback(hObject, eventdata, handles)
global Cloud;

if get(handles.checkbox4,'Value')==1 % execute force GM zero
    Cloud.AngleMagnitude=0;
    set(handles.edit43,'string',num2str(Cloud.AngleMagnitude));% galvo voltage
    if get(handles.checkbox8,'Value') == 1 %yellow stim
        handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
        Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    elseif get(handles.checkbox11,'Value') == 1 %red stim
        handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
        Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    else
        handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
        Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    end
else
    Cloud.AngleMagnitude=Cloud.AngleMagnitudeBackup;
    set(handles.edit43,'string',num2str(Cloud.AngleMagnitude));% galvo voltage
end
end

% --- Executes RUN STIM
function pushbutton2_Callback(hObject, eventdata, handles)
global DAQstate;
global status;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
status = 1;
ToSave.DMDSTATE = 'RUN Stim Test one spot';
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
Stim.Npeaks = str2num(get(handles.edit6,'string')); % Number of blue light pulses in test
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
if get(handles.checkbox2,'Value') == 0
    select = Stim.UT<handles.Setup.Scorepikes.sealtestduration;    
    tempOutput=zeros(size(Stim.Output));
    tempOutput(select,8)=handles.Setup.Scorepikes.sealtestvalue;% seal test voltage:5mV
    queueOutputData(handles.Setup.Daq,tempOutput);
    startForeground(handles.Setup.Daq);
    pause(0.5);
else
    Stim.Output(:,8) =0;
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
% while Stim.Output(find(Stim.UUT>Stim.DelayBorder,1,'first'),1)~=0
%     Stim.DelayBorder=Stim.DelayBorder+0.0001;
% end 
Stim.Array = zeros(size(Stim.UUT));
Stim.Array1=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);%remove last rising edge
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.nonzerovalues = find(Stim.Output(:,2));
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
    Stim.nonzerovalues = find(Stim.Output(:,1));
end
Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
axes(handles.axes3);hold off; cla;
Stim.Output(end,:)=0;

% camera------------------------
if get(handles.checkbox17,'Value') == 1
    handles.Setup.src.TriggerMode = 'Trigger first'; % continue capture
    % handles.Setup.camera.FramesPerTrigger = Inf;
    triggerconfig(handles.Setup.camera, 'hardware', 'Falling edge', 'Extern');
    handles.Setup.src.Exposure = str2double(get(handles.edit58,'string'));
    handles.Setup.src.AutoContrast = 'OFF';
    handles.Setup.src.ExposeOutMode = 'First Row';
    handles.Setup.src.PortSpeedGain = 'Port0-Speed1-100MHz-16bit-Gain1-HDR';
    handles.Setup.src.ClearMode = 'Post-Sequence';
    handles.Setup.camera.FramesPerTrigger = floor(Stim.DelayBorder*1000/handles.Setup.src.Exposure);
    handles.Setup.camera.ROIPosition = [0 0 1200 1200];
    temp=function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),0,Stim.TargetRadius);
    DMDimage=uint8(gather(temp(:,:,1)))*255;
    Stim.Array=circshift(Stim.Array',-round(Stim.DelayBorder*0.99*handles.Setup.Daq.Rate));
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array;
    Stim.VoltageCa = str2num(get(handles.edit63,'string'));
    Stim.Output(:,2) = Stim.VoltageCa;
%     Stim.CropMask2=zeros(size(Stim.UUT))+double(Stim.UUT>Stim.DelayBorder+Stim.DurationMS/1000+(Stim.Npeaks-1)/Stim.FreqHZ);
%     Stim.Output(:,5)=Stim.CropMask2';
    Stim.CropMask2=zeros(size(Stim.UUT))+double(Stim.UUT>Stim.DelayBorder).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(Stim.Npeaks-1)/Stim.FreqHZ);
    Stim.Output(:,5)=circshift(ones(size(Stim.CropMask2'))-Stim.CropMask2',-round(Stim.DelayBorder*0.4*handles.Setup.Daq.Rate)); % trigger camera
    
    if get(handles.checkbox4,'Value') == 1
        Stim.Output(:,3)=handles.Setup.PointCloud.GalvoOffsetVoltage(1);
        Stim.Output(:,4)=handles.Setup.PointCloud.GalvoOffsetVoltage(2);
    else
        Stim.Output(:,3)=Stim.Output(:,3).*Stim.Array'+...
            Stim.CropMask2'.*handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1);
        Stim.Output(:,4)=Stim.Output(:,4).*Stim.Array'+...
            Stim.CropMask2'.*handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2);
    end
    ProjMode = function_CheckProjMode_DMD(handles.Setup);
    if ProjMode==2302
        [handles.Setup]=function_StopProj_DMD(handles.Setup);
        [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
    end
    function_directfeed_DMD(handles.Setup, DMDimage);
    outputSingleScan(handles.Setup.Daq,Stim.Output(1,:));
    pause(30); % wait the yellow laser to be stable;
else
    ProjMode = function_CheckProjMode_DMD(handles.Setup);
    if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
    end
end

for i = 1:floor(str2num(get(handles.edit9,'string')));
    if status == 0; disp('Procedure interrupted');  break; end;
    if get(handles.checkbox17,'Value') == 1
        start(handles.Setup.camera);
    else
        Result.DMDFrames =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),Stim.TargetRadius);
        function_feed_DMD(handles.Setup, Result.DMDFrames);
    end
    queueOutputData(handles.Setup.Daq,Stim.Output);
    Data=startForeground(handles.Setup.Daq);
    if get(handles.checkbox17,'Value') == 1
        stop(handles.Setup.camera);
        currentimage=squeeze(getdata(handles.Setup.camera));
        if i==1
            Images=currentimage;
        else
            Images=cat(3,Images,currentimage);
        end
        mask=zeros(size(currentimage));
        mask(600-Stim.TargetRadius:600+Stim.TargetRadius,600-Stim.TargetRadius:600+Stim.TargetRadius,:)=1;
        currentimagenon0=nonzeros(double(currentimage).*mask);
        DataSave.Image{i}=reshape(currentimagenon0,[numel(currentimagenon0)/size(currentimage,3),size(currentimage,3)]);
        DataSave.Iref{i}=Data;
        Iref = medfilt1(nonzeros(Data.*Stim.Output(:,5)),handles.Setup.Daq.Rate*handles.Setup.src.Exposure/1000);
        if i==1
            Ica=mean(DataSave.Image{i},1);
            Iref1=Iref;
        else
            Ica=cat(2,Ica,mean(DataSave.Image{i},1));
            Iref1=cat(1,Iref1,Iref);
        end
        T_ds=0:handles.Setup.src.Exposure:(numel(Ica)-1)*handles.Setup.src.Exposure;
        T_ref=0:1/handles.Setup.Daq.Rate*1000:(numel(Iref1)-1)/handles.Setup.Daq.Rate*1000;
        axes(handles.axes2);plot(T_ds,Ica);xlabel('Time [ms]'); ylabel('Intensity'); hold off;title('Raw Intensty');
        axes(handles.axes5);plot(T_ref,Iref1);xlabel('Time [ms]'); ylabel('Voltage (V)'); hold off;title('Ref Intensty');
    else
        [ HPData, ~] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        Stim.baseline=mean(Data(1:100));
        [score,~] = function_scorespikes(handles.Setup,Stim, Data);
        axes(handles.axes2);
        plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage'); hold off;title([ 'Ephys Data , Score = ' num2str(floor(score*100)/100)]);
        axes(handles.axes3);
        plot(Stim.UUT,HPData);xlabel('Time s'); ylabel('Voltage'); title([ 'High Passed, Score = ' num2str(floor(score*100)/100)]);
        axes(handles.axes5);
        scatter(i,score,'filled','red'), xlabel('Trial #'); ylabel('score [AU]'); hold on;
        pause(str2num(get(handles.edit10,'string')));
        DataSave{i} = Data;
    end
    axes(handles.axes4);
    if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
        plot(Stim.UUT,Stim.Output(:,2));xlabel('Time [s]'); ylabel('Stim laser [V]');title(['Sequence ' int2str(i) ' of ' int2str(floor(str2num(get(handles.edit9,'string'))))]);
    else
        plot(Stim.UUT,Stim.Output(:,1));xlabel('Time [s]'); ylabel('Stim laser [V]');title(['Sequence ' int2str(i) ' of ' int2str(floor(str2num(get(handles.edit9,'string'))))]);
    end  
    disp(['Repetition #' num2str(i) '/' get(handles.edit9,'string') ' finished!']);
end
if get(handles.checkbox17,'Value') == 1
    CaTrace=Ica-100;%100: background intensity
    dF=zeros(size(CaTrace));
    Stim.cutoffFreq=1/Stim.DelayBorder;
    w=round((1000/handles.Setup.src.Exposure)/Stim.cutoffFreq);
    for ii=1:numel(CaTrace)
        if ii<=w
            F0=min(CaTrace(1:ii+w));
        elseif ii>numel(CaTrace)-w
            F0=min(CaTrace(ii-w:end));
        else
            F0=min(CaTrace(ii-w:ii+w));
        end
        dF(ii)=(CaTrace(ii)-F0)/F0;
    end
    axes(handles.axes3);
    plot(T_ds,dF);xlabel('Time ms'); ylabel('\DeltaF/F'); title('\DeltaF/F');
    figure(401);imshow3D(Images);
    Stim.Camera=handles.Setup.camera;
    Stim.CameraSrc=handles.Setup.src;
end
handles.Setup=function_StopProj_DMD(handles.Setup);
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
ToSave.status = status;
ToSave.type = 'RUN Stim Test';
ToSave.Stim = Stim;
ToSave.Data = DataSave;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave','-v7.3');
    disp('Data Saved')
end
end

% --- Executes Optimize voltage
% function pushbutton3_Callback(hObject, eventdata, handles)
% global DAQstate; global status;status = 1;
% global ToSave;
% global Cloud;
% global SaveID; SaveID=SaveID+1;
% DAQstate = [0 0 0 0 0 0];
% outputSingleScan(handles.Setup.Daq,DAQstate);
% set(handles.slider3,'Value',0)
% set(handles.slider2,'Value',0)
% set(handles.edit3,'string',0);
% set(handles.edit4,'string',0);
% Stim.Npeaks = str2num(get(handles.edit6,'string')); ; % Number of blue light pulses in test
% Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
% Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
% Stim.Voltage = str2num(get(handles.edit11,'string'));
% Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
% Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
% Stim.DutyCycle = str2num(get(handles.edit29,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
% Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
% Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
% [Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
% handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
% Cloud.AnlgeMagnitudeOffset=[0,0];
% [ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
% Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
% if get(handles.checkbox2,'Value') == 0
%     Stim.Output(:,5)=0;
% %     select = Stim.UT<2*handles.Setup.Scorepikes.sealtestduration.* Stim.UT>handles.Setup.Scorepikes.sealtestduration;
% %     Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
%     select = Stim.UT<handles.Setup.Scorepikes.sealtestduration;    
%     tempOutput=Stim.Output;
%     tempOutput(1:end,:)=0;
%     tempOutput(select,5)=1;
%     queueOutputData(handles.Setup.Daq,tempOutput);
%     startForeground(handles.Setup.Daq);
%     pause(0.5);
% else
%     Stim.Output(:,5) =0;
% end
% Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
% Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
% [Stim.LN,Stim.LX] = size(Stim.Output);
% Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
% Stim.Array = Stim.UUT-Stim.UUT;
% Stim.Array1=Stim.Array;
% for i = 1:Stim.Npeaks
%     Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
%     Stim.Array1=Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ-0.0003);
% end
% Stim.Baseline = Stim.Output(:,1);
% if get(handles.checkbox8,'Value') == 1
%     Stim.Output(:,1)=0;
%     Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
% else   
%     Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
%     Stim.Output(:,2)=0;
% end
% Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
% %g = figure(2);
% DAQstate = [DAQstate(1) DAQstate(2) handles.Setup.PointCloud.GalvoOffsetVoltage(1) handles.Setup.PointCloud.GalvoOffsetVoltage(2) 0 0];
% Result.DMDFrames =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
% [handles.Setup,Cloud.targetid] = function_feed_DMD( handles.Setup,squeeze(uint8(gather(Result.DMDFrames(:,:,1)))*255),0);
% axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
% DataSave = {};
% for i = 1:floor(str2num(get(handles.edit9,'string'))) % Do many repetitions
%     for j = 1:numel(Stim.Voltageramp)
%         if status == 0; disp('Procedure interrupted');  break; end;
%         if get(handles.checkbox8,'Value') == 1
%             Stim.Output(:,2) = Stim.Voltageramp(j)*Stim.Baseline.*Stim.Array';
%         else
%             Stim.Output(:,1) = Stim.Voltageramp(j)*Stim.Baseline.*Stim.Array';
%         end
%         queueOutputData(handles.Setup.Daq,Stim.Output);
%         Data=startForeground(handles.Setup.Daq);
%         if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
%         [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
%         Cell.Powerscore{j,i} = score;
%         DataSave{j,i} = Data;
%         axes(handles.axes2); 
%         scatter(Stim.Voltageramp(j),Cell.Powerscore{j,i},'filled','red'); hold on;xlabel('Voltage on laser'); ylabel('Ephys score [AU]');title(['Calibration pass ' int2str(i) 'of' int2str(floor(str2num(get(handles.edit9,'string'))))])
%         axes(handles.axes3);
%         plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage'); title([ 'Score = ' num2str(floor(Cell.Powerscore{j,i}*100)/100)]);
%         axes(handles.axes4);
%         plot(Stim.UUT,Stim.Output(:,1));xlabel('Time [s]'); ylabel('Stim laser [V]');title(['Sequence ' int2str(i) ' of ' int2str(floor(str2num(get(handles.edit9,'string'))))]);
% 
%     end
% end
% if status == 0
% else;
%     Cell.BestLaserVoltage = ginput(1); %close(g);
%     scatter(Cell.BestLaserVoltage(1),Cell.BestLaserVoltage(2),'filled','blue'); title('Optimized voltage selected');
%     Cell.BestLaserVoltage = Cell.BestLaserVoltage(1);
%     set(handles.edit11,'string',num2str(Cell.BestLaserVoltage));
%     ToSave.type = 'Voltage Optimization';
%     ToSave.Stim = Stim;
%     ToSave.Data = DataSave;
%     ToSave.status = status;
%     ToSave.VoltageScores = Cell.Powerscore;
%     ToSave.VoltageRamp = Stim.Voltageramp;
%     if get(handles.checkbox1,'Value') == 1
%         filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
%         save(filename,'ToSave');
%         disp('Data Saved')
%     end  
% end
% end

%--------------optimize voltage------------------
function pushbutton3_Callback(hObject, eventdata, handles)
global DAQstate;
global status;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
status = 1;
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
% Stim.Npeaks = str2num(get(handles.edit6,'string')); % Number of blue light pulses in test
Stim.Npeaks = 1;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.Output(:,5) =0;
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array  = Stim.Array +double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);%remove last rising edge
end
Stim.Baseline = repmat(double(Stim.Subclock<0.8),[1 Stim.NumberCycles])';
if str2num(get(handles.edit18,'string')) == 0 %z=0
    Stim.Output(:,6)=Stim.Array1';
else
    Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
end
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
% Stim.Output(:,6)=Stim.Array1';    
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
Stim.Output(end,:)=0;
DataSave = {};
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end

for i = 1:floor(str2num(get(handles.edit9,'string'))) % Do many repetitions
    for j = 1:numel(Stim.Voltageramp)
        if status == 0; disp('Procedure interrupted');  break; end;
        if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1 
            if str2num(get(handles.edit18,'string')) == 0
                 Stim.Output(:,2) = Stim.Voltageramp(j).*Stim.Array';
            else
                 Stim.Output(:,2) = Stim.Voltageramp(j)*Stim.Baseline.*Stim.Array';
            end
            Stim.Output(:,1)=0;
            Stim.nonzerovalues = find(Stim.Output(:,2));
        else
            if str2num(get(handles.edit18,'string')) == 0
                 Stim.Output(:,1) = Stim.Voltageramp(j).*Stim.Array';
            else
                 Stim.Output(:,1) = Stim.Voltageramp(j)*Stim.Baseline.*Stim.Array';
            end
            Stim.Output(:,2)=0;
            Stim.nonzerovalues = find(Stim.Output(:,1));
        end
        if j==1
            if get(handles.checkbox13,'Value') == 0
                Result.DMDFrames =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
            else
                Stim.Numoffoci=max(str2num(get(handles.edit57,'string')),str2num(get(handles.edit56,'string')));
                X=UX(randi(numel(UX),[Stim.Numoffoci,1]));
                Y=UY(randi(numel(UY),[Stim.Numoffoci,1]));
                Z=UZ(randi(numel(UZ),[Stim.Numoffoci,1]));
                Result.DMDFrames = function_makespots_ori(handles.Setup,X,Y,Z,Stim.TargetRadius*ones(size(X)));
            end
        end
        if str2num(get(handles.edit18,'string')) == 0
            [handles.Setup,~]=function_feed_DMD(handles.Setup, uint8(gather(Result.DMDFrames(:,:,1)))*255);
        else
            [handles.Setup,~]=function_feed_DMD(handles.Setup, uint8(gather(Result.DMDFrames))*255);
        end
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        Stim.baseline=mean(Data(1:100));
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
       
        Cell.Powerscore{j,i} = score;
        DataSave{j,i} = Data;
        axes(handles.axes2); 
        scatter(Stim.Voltageramp(j),Cell.Powerscore{j,i},'filled','red'); hold on;xlabel('Voltage on laser'); ylabel('Ephys score [AU]');title(['Calibration pass ' int2str(i) 'of' int2str(floor(str2num(get(handles.edit9,'string'))))])
        axes(handles.axes3);
        plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage'); title([ 'Score = ' num2str(floor(Cell.Powerscore{j,i}*100)/100)]);
        axes(handles.axes4);
        if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
        plot(Stim.UUT,Stim.Output(:,2));xlabel('Time [s]'); ylabel('Stim laser [V]');title(['Sequence ' int2str(i) ' of ' int2str(floor(str2num(get(handles.edit9,'string'))))]);
        else
          plot(Stim.UUT,Stim.Output(:,1));xlabel('Time [s]'); ylabel('Stim laser [V]');title(['Sequence ' int2str(i) ' of ' int2str(floor(str2num(get(handles.edit9,'string'))))]);
        end  
    end
end
if status == 0
else;
    Cell.BestLaserVoltage = ginput(1); %close(g);
    scatter(Cell.BestLaserVoltage(1),Cell.BestLaserVoltage(2),'filled','blue'); title('Optimized voltage selected');
    Cell.BestLaserVoltage = Cell.BestLaserVoltage(1);
    set(handles.edit13,'string',num2str(Cell.BestLaserVoltage));
    set(handles.edit11,'string',num2str(Cell.BestLaserVoltage));
    ToSave.type = 'Voltage Optimization';
    ToSave.Stim = Stim;
    ToSave.Data = DataSave;
    ToSave.status = status;
    ToSave.VoltageScores = Cell.Powerscore;
    ToSave.VoltageRamp = Stim.Voltageramp;
    if get(handles.checkbox1,'Value') == 1
        filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
        save(filename,'ToSave');
        disp('Data Saved')
    end  
end
end


% --- Executes DARK DMD
function pushbutton5_Callback(hObject, eventdata, handles)
global ToSave;
ToSave.DMDSTATE = 'DMD Dark';
Stim.BlankFrame = zeros(handles.Setup.DMD.LX,handles.Setup.DMD.LY);
ProjMode = function_CheckProjMode_DMD(Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(Setup,'master');
end
function_directfeed_DMD( handles.Setup,Stim.BlankFrame );
global DAQstate;
outputSingleScan(handles.Setup.Daq,DAQstate);
end

% --- Executes "target on DMD"
function pushbutton6_Callback(hObject, eventdata, handles)
global ToSave;
global Cloud;
global DAQstate;
ToSave.DMDSTATE = 'DMD Showing one target';
handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
DAQstate = [DAQstate(1) DAQstate(2) handles.Setup.PointCloud.GalvoOffsetVoltage(1) handles.Setup.PointCloud.GalvoOffsetVoltage(2) 0 0 1 0];
Result.DMDFrames =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
end
[handles.Setup,Cloud.targetid] = function_feed_DMD( handles.Setup,squeeze(uint8(gather(Result.DMDFrames(:,:,1)))*255));
for i=1:4
    DAQstate(6)=rem(i,2);
    outputSingleScan(handles.Setup.Daq,DAQstate);
end
end

% --- Executes BRIGHT DMD
function pushbutton7_Callback(hObject, eventdata, handles)
global ToSave;
global Cloud;
global DAQstate;
ToSave.DMDSTATE = 'DMD Bright';
Stim.BrightFrame = ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY, 'uint8')*255;
handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
DAQstate = [DAQstate(1) DAQstate(2) handles.Setup.PointCloud.GalvoOffsetVoltage(1) handles.Setup.PointCloud.GalvoOffsetVoltage(2) 0 0 1 0];
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
end
[handles.Setup,Cloud.targetid] = function_feed_DMD( handles.Setup,Stim.BrightFrame);
for i=1:4
    DAQstate(6)=rem(i,2);
    outputSingleScan(handles.Setup.Daq,DAQstate);
end
end


% --- Executes PPSF Mechanical
% function pushbutton17_Callback(hObject, eventdata, handles)
% global DAQstate; global status;status = 1;
% global StageStatus;
% global Stage;
% global Cloud;
% if StageStatus==0
%     Stage  = function_Start_stage( handles.Setup.SutterStage );
%     StageStatus = 1;
% end
% function_Zero_stage( Stage );
% DAQstate = [0 0 0 0 0 0];
% global ToSave;
% global SaveID; SaveID=SaveID+1;
% ToSave.DMDSTATE = 'DMD Digital PPSF Sequence';
% outputSingleScan(handles.Setup.Daq,DAQstate);
% set(handles.slider3,'Value',0);
% set(handles.slider2,'Value',0);
% set(handles.edit3,'string',0);
% set(handles.edit4,'string',0);
% Stim.Npeaks = floor(str2num(get(handles.edit6,'string'))); % Number of blue light pulses in test
% Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
% Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
% Stim.Voltage = str2num(get(handles.edit11,'string'));
% Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
% Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
% Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
% Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
% Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
% [Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
% if get(handles.checkbox8,'Value') == 0
%     handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
%     Cloud.AnlgeMagnitudeOffset=[0,0];
% else
%     handles.Setup.PointCloud.GalvoOffsetVoltage=[0.2,2.7];
%     Cloud.AnlgeMagnitudeOffset=[0.2,2.7];
% end
% [ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
% Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
% if get(handles.checkbox2,'Value') == 0
%     Stim.Output(:,5)=0;
% %     select = Stim.UT<2*handles.Setup.Scorepikes.sealtestduration.* Stim.UT>handles.Setup.Scorepikes.sealtestduration;
% %     Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
%     select = Stim.UT<handles.Setup.Scorepikes.sealtestduration;    
%     tempOutput=Stim.Output;
%     tempOutput(1:end,:)=0;
%     tempOutput(select,5)=1;
%     queueOutputData(handles.Setup.Daq,tempOutput);
%     startForeground(handles.Setup.Daq);
%     pause(0.5);
% else
%     Stim.Output(:,5) =0;
% end
% Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
% Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
% [Stim.LN,Stim.LX] = size(Stim.Output);
% Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
% Stim.Array = Stim.UUT-Stim.UUT;
% for i = 1:Stim.Npeaks
%     Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
% end
% Stim.Baseline = Stim.Output(:,1);
% if get(handles.checkbox8,'Value') == 1
%     Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
%     Stim.nonzerovalues = find(Stim.Output(:,2));
%     Stim.Output(:,1)=0;
% else   
%     Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
%     Stim.Output(:,2)=0;
%     Stim.nonzerovalues = find(Stim.Output(:,1));
% end
% 
% if get(handles.checkbox3,'Value') == 1
% UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
% UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
% UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
% else
% UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
% UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
% UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
% end
% DataSave.UX = UX;
% DataSave.UY = UY;
% DataSave.UZ = UZ;
% DataSave.X = {};  DataSave.Y = {}; DataSave.Z = {}; DataSave.SX = {}; DataSave.SY = {}; DataSave.SZ = {};
% Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
% Cloud.divider = floor(str2num(get(handles.edit42,'string')));
% Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
% Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
% Cloud.divider = floor(str2num(get(handles.edit42,'string')));
% Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
% handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
% Result.DMDFrames =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
% function_feed_DMD(uint8(gather(Result.DMDFrames))*255);
% axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
% Stim.Output(end,:)=0;
% for j = 1:floor(str2num(get(handles.edit9,'string'))); %do as many repetitions as needed
%     if status == 0; disp('Procedure interrupted'); break; end;
%     for i = 1:floor(str2num(get(handles.edit23,'string')));
%         if status == 0; disp('Procedure interrupted'); break; end;
%         Position = [UX(i) 0 0];
%         function_Goto_stage( Stage,Position );
%         set(handles.edit26,'string',num2str(Position(1)));%Stage X
%         set(handles.edit27,'string',num2str(Position(2)));%Stage Y
%         set(handles.edit28,'string',num2str(Position(3)));pause(1);%Stage Z
%         queueOutputData(handles.Setup.Daq,Stim.Output);
%         Data=startForeground(handles.Setup.Daq);
%         [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
%         if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
%         [score, odata] = function_scorespikes(handles.Setup,Stim, Data );
%         axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys data , Score = ' num2str(score)]);
%         axes(handles.axes2);scatter(UX(i), score,'red','filled'); hold on; xlabel('X pixels'); ylabel('Score'); title('X PPSF');
%         pause(str2num(get(handles.edit10,'string')));
%         DataSave.X{i,j} = Data;
%         DataSave.SX{i,j} = score;
%     end
%     for i = 1:floor(str2num(get(handles.edit24,'string')));
%         if status == 0; disp('Procedure interrupted'); break; end;
%         Position = [0 UY(i) 0];
%         function_Goto_stage( Stage,Position );
%         set(handles.edit26,'string',num2str(Position(1)));%Stage X
%         set(handles.edit27,'string',num2str(Position(2)));%Stage Y
%         set(handles.edit28,'string',num2str(Position(3)));pause(1);%Stage Z
%         queueOutputData(handles.Setup.Daq,Stim.Output);
%         Data=startForeground(handles.Setup.Daq);
%         [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
%         if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
%         [score, odata] = function_scorespikes(handles.Setup,Stim, Data );
%         axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys data , Score = ' num2str(score)]);
%         axes(handles.axes3);scatter(UY(i), score,'red','filled'); hold on;xlabel('Y pixels'); ylabel('Score'); title('Y PPSF');
%         pause(str2num(get(handles.edit10,'string')));
%         DataSave.Y{i,j} = Data;
%         DataSave.SY{i,j} = score;
%     end
%     for i = 1:floor(str2num(get(handles.edit25,'string')));
%         if status == 0;disp('Procedure interrupted');break; end;
%         Position = [ 0 0 UZ(i)];
%         function_Goto_stage( Stage,Position );
%         set(handles.edit26,'string',num2str(Position(1)));%Stage X
%         set(handles.edit27,'string',num2str(Position(2)));%Stage Y
%         set(handles.edit28,'string',num2str(Position(3)));pause(1);%Stage Z
%         queueOutputData(handles.Setup.Daq,Stim.Output);
%         Data=startForeground(handles.Setup.Daq);
%         [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
%         if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
%         [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
%         axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys Data , Score = ' num2str(score)]);
%         axes(handles.axes4); scatter(UZ(i), score,'red','filled'); hold on;xlabel('Z pixels'); ylabel('Score'); title('Z PPSF');
%         pause(str2num(get(handles.edit10,'string')));
%         DataSave.Z{i,j} = Data;
%         DataSave.SZ{i,j} = score;
%     end
% end
% Position = [0 0 0];
% function_Goto_stage( Stage,Position );
% set(handles.edit26,'string',num2str(Position(1)));%Stage X
% set(handles.edit27,'string',num2str(Position(2)));%Stage Y
% set(handles.edit28,'string',num2str(Position(3)));pause(1);%Stage Z
% ToSave.type = 'Mechanical PPSF';
% ToSave.Stim = Stim;
% ToSave.status = status;
% ToSave.Data = DataSave;
% if get(handles.checkbox1,'Value') == 1
%     filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
%     save(filename,'ToSave');
%     disp('Data Saved')
% end
% end


% --- Executes Digital PPSF with just full XZ plane
function pushbutton22_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'DMD Digital PPSF Sequence';
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
Stim.Npeaks = floor(str2num(get(handles.edit6,'string'))); % Number of blue light pulses in test
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
if get(handles.checkbox2,'Value') == 0
    Stim.Output(:,5)=0;
%     select = Stim.UT<2*handles.Setup.Scorepikes.sealtestduration.* Stim.UT>handles.Setup.Scorepikes.sealtestduration;
%     Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
    select = Stim.UT<handles.Setup.Scorepikes.sealtestduration;    
    tempOutput=Stim.Output;
    tempOutput(1:end,:)=0;
    tempOutput(select,5)=1;
    queueOutputData(handles.Setup.Daq,tempOutput);
    startForeground(handles.Setup.Daq);
    pause(0.1);
else
    Stim.Output(:,5) =0;
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.nonzerovalues = find(Stim.Output(:,2));
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
    Stim.nonzerovalues = find(Stim.Output(:,1));
end
Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UZ = UZ;
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end
[PXX,PZZ] = ndgrid(UX,UZ);
DataSave.PXX = PXX(:);
DataSave.PZZ = PZZ(:);

DataSave.XZ = {};  DataSave.SXZ = {};
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
for j = 1:floor(str2num(get(handles.edit9,'string'))) %do as many repetitions as needed
    if j==1
        for i=1:numel(PXX(:))
            if status == 0; disp('Procedure interrupted'); break; end;
            [DMDFrames] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string'))+DataSave.PXX(i),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string'))+DataSave.PZZ(i),str2num(get(handles.edit19,'string')));
            [handles.Setup,Cloud.sequenceid2D(i)] = function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames))*255);
        end
    end
    handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
    for i=1:numel(PXX(:))
        if status == 0; disp('Procedure interrupted'); break; end;
        Currentsequenceid=Cloud.sequenceid2D(i);
        handles.Setup = function_StartProj_DMD(handles.Setup, Currentsequenceid);
%         StimPower=Stim.Voltage+sqrt(PXX(i)^2+PZZ(i)^2)/sqrt(max(UX)^2+max(UZ)^2)*(5-Stim.Voltage);%increase power with distance to center
%         if get(handles.checkbox8,'Value') == 1
%             Stim.Output(:,2)=StimPower*Stim.Baseline.*Stim.Array';
%         else
%             Stim.Output(:,1)=StimPower*Stim.Baseline.*Stim.Array';
%         end
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
        Stim.baseline=mean(Data(1:100));
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else handles.Setup.Scorepikes.Method=0; end;
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
        axes(handles.axes2);  scatter(DataSave.PXX(i),DataSave.PZZ(i), 25,score,'filled'); colorbar; hold on; xlabel('X pixels'); ylabel('Z pixels'); title('XZ PPSF score'); 
        axes(handles.axes5);  plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage');
        DataSave.XZ{i,j} = Data;
        DataSave.SXZ{i,j} = score;
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
    end
    disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
end
try
ToSave.type = 'Digital PPSF Full XZ plane';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
catch; end;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end
end


function edit42_Callback(hObject, eventdata, handles)
global Cloud;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
end

function edit43_Callback(hObject, eventdata, handles)
global Cloud;
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
if Cloud.AngleMagnitude>0
Cloud.AngleMagnitudeBackup = Cloud.AngleMagnitude;
end
end

function edit44_Callback(hObject, eventdata, handles)
global Cloud;
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
handles.Setup.PointCloud.CycleLength=Cloud.CycleLength;
end


% --- Executes Make PPSF digitally
function pushbutton8_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'DMD Digital PPSF Sequence';
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
Stim.Npeaks = floor(str2num(get(handles.edit6,'string'))); % Number of blue light pulses in test
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
if get(handles.checkbox2,'Value') == 0
    Stim.Output(:,5)=0;
%     select = Stim.UT<2*handles.Setup.Scorepikes.sealtestduration.* Stim.UT>handles.Setup.Scorepikes.sealtestduration;
%     Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
    select = Stim.UT<handles.Setup.Scorepikes.sealtestduration;    
    tempOutput=Stim.Output;
    tempOutput(1:end,:)=0;
    tempOutput(select,5)=1;
    queueOutputData(handles.Setup.Daq,tempOutput);
    startForeground(handles.Setup.Daq);
    pause(0.5);
else
    Stim.Output(:,5) =0;
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
Stim.Array2=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array2 = Stim.Array2+double(Stim.UUT>Stim.DelayBorder*0.98+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.nonzerovalues = find(Stim.Output(:,2));
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
    Stim.nonzerovalues = find(Stim.Output(:,1));
end
Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
DataSave.X = {};  DataSave.Y = {}; DataSave.Z = {}; 
DataSave.SX = {}; DataSave.SY = {}; DataSave.SZ = {};
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
Stim.Output(end,:)=0;
Nx=floor(str2num(get(handles.edit23,'string')));
Ny=floor(str2num(get(handles.edit24,'string')));
Nz=floor(str2num(get(handles.edit25,'string')));
Cloud.sequenceid1D.x=zeros(Nx,1);
Cloud.sequenceid1D.y=zeros(Ny,1);
Cloud.sequenceid1D.z=zeros(Nz,1);
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

handles.Setup.DMD.SequenceControl.RepeatModeValue=Stim.Npeaks;
for j = 1:floor(str2num(get(handles.edit9,'string'))); %do as many repetitions as needed
    if status == 0; disp('Procedure interrupted'); break; end;
% Generate Patterns
    %generate patterns for x direction 
    if j==1 % only need to generate pattern once
        for i=1:Nx
            if status == 0; disp('Procedure interrupted'); break; end;
            [DMDFrames] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string'))+UX(i),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
            [handles.Setup,Cloud.sequenceid1D.x(i)] = function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames))*255);
            
        end
        for i=1:Ny
            if status == 0; disp('Procedure interrupted'); break; end;
            [DMDFrames] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string'))+UY(i),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
            if get(handles.checkbox17,'Value') == 1
                temp=function_makespots(handles.Setup,0,0,0,str2num(get(handles.edit19,'string')));
                DMDFrames = cat(3,DMDFrames, temp(:,:,1));
            end
            [handles.Setup,Cloud.sequenceid1D.y(i)] = function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames))*255);
        end
         for i=1:Nz
            if status == 0; disp('Procedure interrupted'); break; end;
            [DMDFrames] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string'))+UZ(i),str2num(get(handles.edit19,'string')));
            if get(handles.checkbox17,'Value') == 1
                temp=function_makespots(handles.Setup,0,0,0,str2num(get(handles.edit19,'string')));
                DMDFrames = cat(3,DMDFrames, temp(:,:,1));
            end
            [handles.Setup,Cloud.sequenceid1D.z(i)] = function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames))*255);
         end
    end
        
    for i=1:Nx
        if status == 0; disp('Procedure interrupted'); break; end;
        if get(handles.checkbox17,'Value') == 0
            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.sequenceid1D.x(i));
            queueOutputData(handles.Setup.Daq,Stim.Output);
            Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
            if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else handles.Setup.Scorepikes.Method=0; end;
            Stim.baseline=mean(Data(end-100:end));
            [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
            axes(handles.axes2); scatter(UX(i), score,'red','filled'); hold on; xlabel('X pixels'); ylabel('Score'); title('X PPSF');
            axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage');
            DataSave.X{i,j} = Data;
            DataSave.SX{i,j} = score;
        else
            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.sequenceid1D.x(i));
%             queueOutputData(handles.Setup.Daq,Stim.Output);
            queueOutputData(handles.Setup.Daq,Stim.Output(1:round(Stim.Npeaks/Stim.FreqHZ*handles.Setup.Daq.Rate),:));
            startForeground(handles.Setup.Daq);
            
            [handles.Setup,~] = function_feed_DMD( handles.Setup,ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY, 'uint8')*255);
            start(handles.Setup.camera);
            queueOutputData(handles.Setup.Daq,Stim.Output(round(Stim.Npeaks/Stim.FreqHZ*handles.Setup.Daq.Rate)+1:end,:));
            Iref=startForeground(handles.Setup.Daq);
            Images=squeeze(getdata(handles.Setup.camera));
            stop(handles.Setup.camera);

            % plot the result
            Mask = ones(Stim.ROIoffset(3),Stim.ROIoffset(4));
            Ica = function_IntensityTrace(Images,Mask);
            score=function_scoreIntensityTrace(Stim,Ica);
            Ireffilter=medfilt1(Iref,handles.Setup.Daq.Rate/100);
            axes(handles.axes2); scatter(UX(i), score,'red','filled'); hold on; xlabel('X pixels'); ylabel('Score (Intensity)'); title('X PPSF');
            axes(handles.axes5);plot(Ica);title('flurescence of ROIs');xlabel('frames');
            axes(handles.axes3);imagesc(mean(Images,3));colorbar;title('image of ROIs')
%             axes(handles.axes4);plot(Stim.UUT, Ireffilter);title('Ref: laser intensity');xlabel('time (s)');
            DataSave.X{i,j} = Images;
            DataSave.SX{i,j} = score;
%             DataSave.Iref{i,j} = Iref;
        end
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
    end
    pause(str2num(get(handles.edit10,'string')));
    
    for i=1:Ny
        if get(handles.checkbox17,'Value') == 0
            if status == 0; disp('Procedure interrupted'); break; end;
            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.sequenceid1D.y(i));
            queueOutputData(handles.Setup.Daq,Stim.Output);
            Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
            if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else handles.Setup.Scorepikes.Method=0; end;
            Stim.baseline=mean(Data(1:100));
            [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
            axes(handles.axes2); scatter(UY(i), score,'red','filled'); hold on; xlabel('Y pixels'); ylabel('Score'); title('Y PPSF');
            axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage');
            DataSave.Y{i,j} = Data;
            DataSave.SY{i,j} = score;
        else
            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.sequenceid1D.y(i));
            start(handles.Setup.camera);
            queueOutputData(handles.Setup.Daq,Stim.Output);
            Iref=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
            Images=squeeze(getdata(handles.Setup.camera));
%             if j==1
%                 Images = Imagetemp;
%             else
%                 Images = cat(3, Images, Imagetemp);
%             end
            stop(handles.Setup.camera);

            % plot the result
            Mask = ones(Stim.ROIoffset(3),Stim.ROIoffset(4));
            Ica = function_IntensityTrace(Images,Mask);
            score=function_scoreIntensityTrace(Stim,Ica);
           
            axes(handles.axes2); scatter(UY(i), score,'red','filled'); hold on; xlabel('Y pixels'); ylabel('Score (Intensity)'); title('Y PPSF');
            axes(handles.axes5);plot(Ica);title('flurescence of ROIs');xlabel('frames');
            axes(handles.axes3);imagesc(mean(Images,3));colorbar;title('image of ROIs')
            axes(handles.axes4);plot(Stim.UUT, medfilt1(Iref,handles.Setup.Daq.Rate/80));title('Ref: laser intensity');xlabel('time (s)');
            DataSave.Y{i,j} = Images;
            DataSave.SY{i,j} = score;
            DataSave.IrefY{i,j} = Iref;
        end
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
        
    end
        pause(str2num(get(handles.edit10,'string')));
        
    for i=1:Nz
        if get(handles.checkbox17,'Value') == 0
            if status == 0; disp('Procedure interrupted'); break; end;
            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.sequenceid1D.z(i));
            queueOutputData(handles.Setup.Daq,Stim.Output);
            Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
            if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else handles.Setup.Scorepikes.Method=0; end;
            Stim.baseline=mean(Data(1:100));
            [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
            axes(handles.axes3); scatter(UZ(i), score,'red','filled'); hold on; xlabel('Z pixels'); ylabel('Score'); title('Z PPSF');
            axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage');
            DataSave.Z{i,j} = Data;        DataSave.SZ{i,j} = score;
            
        else
            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.sequenceid1D.z(i));
            start(handles.Setup.camera);
            queueOutputData(handles.Setup.Daq,Stim.Output);
            Iref=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
            Images=squeeze(getdata(handles.Setup.camera));
%             if j==1
%                 Images = Imagetemp;
%             else
%                 Images = cat(3, Images, Imagetemp);
%             end
            stop(handles.Setup.camera);

            % plot the result
            Mask = ones(Stim.ROIoffset(3),Stim.ROIoffset(4));
            Ica = function_IntensityTrace(Images,Mask);
            score=function_scoreIntensityTrace(Stim,Ica);
            axes(handles.axes2); scatter(UZ(i), score,'red','filled'); hold on; xlabel('Z pixels'); ylabel('Score (Intensity)'); title('Z PPSF');
            axes(handles.axes5);plot(Ica);title('flurescence of ROIs');xlabel('frames');
            axes(handles.axes3);imagesc(mean(Images,3));colorbar;title('image of ROIs')
            axes(handles.axes4);plot(Stim.UUT, medfilt1(Iref,handles.Setup.Daq.Rate/80));title('Ref: laser intensity');xlabel('time (s)');
            DataSave.Z{i,j} = Images;
            DataSave.SZ{i,j} = score;
            DataSave.IrefZ{i,j} = Iref;
        end
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
    end
    pause(str2num(get(handles.edit10,'string')));
    
    figure(300);set(gcf,'position',[100,100,250,600]);
    subplot(3, 1, 1);errorbar(UX, mean(cell2mat(DataSave.SX),2),std(cell2mat(DataSave.SX),[],2));
    subplot(3, 1, 2);errorbar(UY, mean(cell2mat(DataSave.SY),2),std(cell2mat(DataSave.SY),[],2));
    subplot(3, 1, 3);errorbar(UZ, mean(cell2mat(DataSave.SZ),2),std(cell2mat(DataSave.SZ),[],2));
    disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
end
try
    Stim.Camera=handles.Setup.camera;
    Stim.CameraSrc=handles.Setup.src;
    ToSave.type = 'Digital PPSF';
    ToSave.Stim = Stim;
    ToSave.status = status;
    ToSave.Data = DataSave;
catch;end;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave','-v7.3');
    disp('Data Saved')
end
end

% ---  Executes  Stage on
function pushbutton9_Callback(hObject, eventdata, handles)
global Stage;
global StageStatus;
try
    Stage  = function_Start_stage( handles.Setup.MechStageComport );
    function_Zero_stage( Stage );
    StageStatus = 1;
    set(handles.edit26,'string',num2str(0));%Stage X
    set(handles.edit27,'string',num2str(0));%Stage Y
    set(handles.edit28,'string',num2str(0));%Stage Z
catch
    set(handles.edit26,'string','ERR');%Stage X
    set(handles.edit27,'string','ERR');%Stage Y
    set(handles.edit28,'string','ERR');%Stage Z
end
end

% --- Executes Stage Off
function pushbutton10_Callback(hObject, eventdata, handles)
global Stage;
global StageStatus;
try
    StageStatus=0;
    function_Stop_stage( Stage );
    set(handles.edit26,'string','OFF');%Stage X
    set(handles.edit27,'string','OFF');%Stage Y
    set(handles.edit28,'string','OFF');%Stage Z
catch
    set(handles.edit26,'string','ERR');%Stage X
    set(handles.edit27,'string','ERR');%Stage Y
    set(handles.edit28,'string','ERR');%Stage Z
end

end

% --- Executes Define stage zero
function pushbutton11_Callback(hObject, eventdata, handles)
global Stage;
global StageStatus;
try
    function_Zero_stage( Stage );
    handles.stageposition = [0 0 0];
    set(handles.edit26,'string',num2str(0));%Stage X
    set(handles.edit27,'string',num2str(0));%Stage Y
    set(handles.edit28,'string',num2str(0));%Stage Z
catch
    set(handles.edit26,'string','ERR');%Stage X
    set(handles.edit27,'string','ERR');%Stage Y
    set(handles.edit28,'string','ERR');%Stage Z
end
end

% --- Executes read stage position
function pushbutton12_Callback(hObject, eventdata, handles)
try
    global Stage;
    [ positionmicron ] = function_Get_stage(Stage );
    handles.stageposition = positionmicron;
    set(handles.edit26,'string',num2str(positionmicron(1)));%Stage X
    set(handles.edit27,'string',num2str(positionmicron(2)));%Stage Y
    set(handles.edit28,'string',num2str(positionmicron(3)));%Stage Z
catch
    set(handles.edit26,'string','ERR');%Stage X
    set(handles.edit27,'string','ERR');%Stage Y
    set(handles.edit28,'string','ERR');%Stage Z
end
end

% --- Executes go home as defined by zero
function pushbutton13_Callback(hObject, eventdata, handles)
try
    global Stage;
    position = [0 0 0]; function_Goto_stage( Stage,position ); pause(0.5);
catch
    set(handles.edit26,'string','ERR');%Stage X
    set(handles.edit27,'string','ERR');%Stage Y
    set(handles.edit28,'string','ERR');%Stage Z
end
end

% --- Executes Go to position
function pushbutton14_Callback(hObject, eventdata, handles)
try
    global Stage;
    stageposition = [0 0 0];
    stageposition(1) = str2num(get(handles.edit26,'string'));
    stageposition(2) = str2num(get(handles.edit27,'string'));
    stageposition(3) = str2num(get(handles.edit28,'string'));
    function_Goto_stage( Stage,stageposition ); pause(0.5);
catch
    set(handles.edit26,'string','ERR');%Stage X
    set(handles.edit27,'string','ERR');%Stage Y
    set(handles.edit28,'string','ERR');%Stage Z
end
end

% --- Executes  Properly QUIT GUI
function pushbutton15_Callback(hObject, eventdata, handles)
try
    global Stage;
    function_Stop_stage( Stage );
catch
end
close;
end

% --- Executes Interrupts execution
function Interrupt_Callback(hObject, eventdata, handles)
global status;
status = 0; % 1 if ok to continue
end

% --- Executes Adjust voltage for laser 1
function slider2_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
set(handles.edit3,'string',round(10*h)/10);
global DAQstate;
DAQstate(1) = h;
outputSingleScan(handles.Setup.Daq,DAQstate);
end

% --- Executes  Adjust voltage for Laser 2
function slider3_Callback(hObject, eventdata, handles)
h = get(hObject,'Value');
set(handles.edit4,'string',round(10*h)/10);
global DAQstate;
DAQstate(2) = h;
outputSingleScan(handles.Setup.Daq,DAQstate);
end

% --- Executes All lasers off
function pushbutton1_Callback(hObject, eventdata, handles)
global DAQstate;
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0)
set(handles.slider2,'Value',0)
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
end

%%%%%%%%% UNUSED COMMANDS BELOW
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

function varargout = Controller_v0_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
end

function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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

function edit20_Callback(hObject, eventdata, handles)
end

function edit20_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit21_Callback(hObject, eventdata, handles)
end

function edit21_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit22_Callback(hObject, eventdata, handles)
end

function edit22_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit23_Callback(hObject, eventdata, handles)
end

function edit23_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit24_Callback(hObject, eventdata, handles)
end

function edit24_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit25_Callback(hObject, eventdata, handles)
end

function edit25_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit26_Callback(hObject, eventdata, handles)
end

function edit26_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit27_Callback(hObject, eventdata, handles)
end

function edit27_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit28_Callback(hObject, eventdata, handles)
end

function edit28_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit29_Callback(hObject, eventdata, handles)
end

function edit29_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit30_Callback(hObject, eventdata, handles)
end

function edit30_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function checkbox1_Callback(hObject, eventdata, handles)
end

function edit31_Callback(hObject, eventdata, handles)
end

function edit31_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function checkbox2_Callback(hObject, eventdata, handles)
% check to count spike, uncheck for photocurrent
end

function checkbox3_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end


function edit32_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit33_Callback(hObject, eventdata, handles)
end

function edit33_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit34_Callback(hObject, eventdata, handles)
end

function edit34_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit35_Callback(hObject, eventdata, handles)
end
    function edit35_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    end

function uipanel4_CreateFcn(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton21.
function pushbutton20_Callback(hObject, eventdata, handles)
end

function edit39_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit40_Callback(hObject, eventdata, handles)
end


function edit40_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function edit41_Callback(hObject, eventdata, handles)
end
    
function edit41_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit36_Callback(hObject, eventdata, handles)
end

function edit36_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit37_Callback(hObject, eventdata, handles)
end

function edit37_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit38_Callback(hObject, eventdata, handles)
end
    function edit38_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
    end
function edit42_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function edit43_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit44_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in DigitalPPSFXZ.
function DigitalPPSFXZ_Callback(hObject, eventdata, handles)
% hObject    handle to DigitalPPSFXZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global SaveID; SaveID=SaveID+1;
Xratio=handles.Setup.DataAnalysis.Xratio;
Zratio=handles.Setup.DataAnalysis.Zratio;
UX=ToSave.Data.UX*Xratio;
UZ=ToSave.Data.UZ*Zratio;
for j=1:size(ToSave.Data.X,2)
    for i=1:size(ToSave.Data.X,1)
        [ SX(i,j), ~] = function_scorespikes_averageAllpulses( handles.Setup,ToSave.Stim,cell2mat(ToSave.Data.X(i,j)));
        [ SZ(i,j), ~] = function_scorespikes_averageAllpulses( handles.Setup,ToSave.Stim,cell2mat(ToSave.Data.Z(i,j)));
    end
end
[GaussianX, ~] = createFit(UX', mean(SX,2));
[GaussianZ, ~] = createFit(UZ', mean(SZ,2));
FWHMx=function_FWHMofGaussian( GaussianX, UX );
FWHMz=function_FWHMofGaussian( GaussianZ, UZ );
DataSave.FWHMx_XZ=FWHMx;
DataSave.FWHMz_XZ=FWHMz;
DataSave.GaussianX_XZ=GaussianX;
DataSave.GaussianZ_XZ=GaussianZ;
ToSave.DataAnalysis = DataSave;
xzppsf=figure;
figure(xzppsf);
set(xzppsf, 'Position',  [300, 300, 1000, 400])
subplot(1,2,1);plot(GaussianX, UX, mean(SX,2));title(['X PPSF']);
grid on;xlim([min(UX) max(UX)]);xlabel('X (\mum)');
dim = [0.14 0.6 0.3 0.3];
str = {['FWHMx=' num2str(FWHMx) '\mum']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
subplot(1,2,2);plot(GaussianZ, UZ, mean(SZ,2));title(['Z PPSF']);
grid on;xlim([min(UZ) max(UZ)]);xlabel('Z (\mum)');
dim = [0.58 0.6 0.3 0.3];
str = {['FWHMz=' num2str(FWHMz) '\mum']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
    if get(handles.checkbox1,'Value') == 1
        filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
        save(filename,'ToSave');
        saveas(xzppsf,[filename '_XZPPSF_plot.tif']);
        disp('Data & Plot Saved');
    end
end

% --- Executes on button press in DigitalPPSFFullXZ.
function DigitalPPSFFullXZ_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global SaveID; SaveID=SaveID+1;
spikeremoveflag=0;
Xratio=handles.Setup.DataAnalysis.Xratio;
Zratio=handles.Setup.DataAnalysis.Zratio;
UX=ToSave.Data.UX*Xratio;
UZ=ToSave.Data.UZ*Zratio;
    for j=1:size(ToSave.Data.XZ,2)
        if find(cellfun(@isempty,ToSave.Data.XZ(:,j)))
            break;
        else
            for i=1:size(ToSave.Data.XZ,1)
                [ SXZ(i,j), odata, spikeremoveflag] = function_score_removeSpike( handles.Setup,ToSave.Stim,cell2mat(ToSave.Data.XZ(i,j)),spikeremoveflag);            
            end
        end
    end
SXZ2D=reshape(mean(SXZ,2),[length(UX) length(UZ)]);
[~,Ind]=max(mean(SXZ,2));
if numel(Ind)>1
    [~,ind]=min(abs(Ind-length(SXZ)/2));
    Ind=Ind(ind);
end
[I,J]=ind2sub([length(UX) length(UZ)],Ind);
x=(I-J+1):size(SXZ2D,1);
x(x<1)=[];%delete non-positive index
y=1:length(x);
for k=1:length(x)
    diagSXZ(k)=SXZ2D(x(k),y(k));
end
d=sqrt((UX(1)-UX(2))^2+(UZ(1)-UZ(2))^2);
UXZ=-(J-1)*d:d:(length(x)-J)*d;    
[GaussianXZ, ~] = createFit(UXZ, diagSXZ);
FWHMxz=function_FWHMofGaussian( GaussianXZ, UXZ );
DataSave.FWHMxzDiagonal=FWHMxz;
DataSave.GaussianXZ_XZdiagonal=GaussianXZ;
ToSave.DataAnalysis = DataSave;
xzdiag=figure;
figure(xzdiag);
set(xzdiag, 'Position',  [300, 300,1000, 400])
subplot(1,2,1);
imagesc(UX,UZ,SXZ2D');
title('XZ PPSF');
axis image;xlabel('X (\mum)');ylabel('Z (\mum)');
colorbar;
subplot(1,2,2);plot(GaussianXZ, UXZ, diagSXZ);title('diag XZ PPSF');
grid on;xlim([min(UZ) max(UZ)]);xlabel('XZ (\mum)');
dim = [0.58 0.6 0.3 0.3];
str = {['FWHMxz=' num2str(FWHMxz) '\mum']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
dim = [0.4 0.7 0.3 0.3];
str = {['Repeat ' num2str(size(SXZ,2)) ', spikeremove=' num2str(spikeremoveflag)]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
if get(handles.checkbox1,'Value') == 1
   filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
   save(filename,'ToSave');
   saveas(xzdiag,[filename '_XZdiagonalPPSF_plot.tif']);
   disp('Data & Plot Saved');
end
end

% --- Executes on button press in SpikeCountingPPSF.
function SpikeCountingPPSF_Callback(hObject, eventdata, handles)
% hObject    handle to SpikeCountingPPSF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global SaveID; SaveID=SaveID+1;
Xratio=handles.Setup.DataAnalysis.Xratio;
Zratio=handles.Setup.DataAnalysis.Zratio;
UX=ToSave.Data.UX*Xratio;
UZ=ToSave.Data.UZ*Zratio;
for j=1:size(ToSave.Data.X,2)
        if isempty(cell2mat(ToSave.Data.X(:,j)))
            break;
        else
            for i=1:size(ToSave.Data.X,1)
                [ SX(i,j), ~] = function_scorespikes( handles.Setup,ToSave.Stim,cell2mat(ToSave.Data.X(i,j)));
                [ SZ(i,j), ~] = function_scorespikes( handles.Setup,ToSave.Stim,cell2mat(ToSave.Data.Z(i,j)));
            end
        end
end
[GaussianX, ~] = createFit(UX', mean(SX,2));
[GaussianZ, ~] = createFit(UZ', mean(SZ,2));
FWHMx=function_FWHMofGaussian( GaussianX, UX );
FWHMz=function_FWHMofGaussian( GaussianZ, UZ );
DataSave.FWHMx_Spike=FWHMx;
DataSave.FWHMz_Spike=FWHMz;
DataSave.GaussianX_Spike=GaussianX;
DataSave.GaussianZ_Spike=GaussianZ;
ToSave.DataAnalysis = DataSave;
spikeCount=figure;
figure(spikeCount);
set(spikeCount, 'Position',  [300, 300, 1000, 400])
subplot(1,2,1);
plot(GaussianX, UX, mean(SX,2));
title('X spike PPSF');
grid on;xlim([min(UX) max(UX)]);xlabel('X (\mum)');ylabel('Num of Spike');
dim = [0.14 0.6 0.3 0.3];
str = {['FWHMx=' num2str(FWHMx) '\mum']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
subplot(1,2,2);
plot(GaussianZ, UZ, mean(SZ,2));
title('Z spike PPSF');
grid on;xlim([min(UZ) max(UZ)]);xlabel('Z (\mum)');ylabel('Num of Spike');
dim = [0.58 0.6 0.3 0.3];
str = {['FWHMz=' num2str(FWHMz) '\mum']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
figure(spikeCount);
dim = [0.47 0.68 0.3 0.3];
str = {['Repeat ' num2str(size(ToSave.Data.X,2))]};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
if get(handles.checkbox1,'Value') == 1
   filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
   save(filename,'ToSave');
   saveas(spikeCount,[filename '_SpikeCountPPSF_plot.tif']);
   disp('Data & Plot Saved');
end
end


% --- red stimuation laser.
function checkbox8_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox8
global Cloud;
handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;
end

% % --- 3D PPSF XYZ with preload mask of 15um aperture
% function pushbutton32_Callback(hObject, eventdata, handles)
% global DAQstate; global status;status = 1;
% DAQstate = [0 0 0 0 0 0];
% global ToSave;
% global Cloud;
% global SaveID; SaveID=SaveID+1;
% ToSave.DMDSTATE = 'DMD Digital PPSF Sequence';
% outputSingleScan(handles.Setup.Daq,DAQstate);
% set(handles.slider3,'Value',0);
% set(handles.slider2,'Value',0);
% set(handles.edit3,'string',0);
% set(handles.edit4,'string',0);
% 
% UX = linspace(-100,100,11);%range should be changed
% UY = linspace(-100,100,11);
% UZ = linspace(-100,100,21);
% 
% DataSave.UX = UX;
% DataSave.UY = UY;
% DataSave.UZ = UZ;
% Nx=11;
% Ny=11;
% Nz=21;
% 
% Stim.Npeaks = 1; %force speak to be 1
% Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
% Stim.FreqHZ= 1000/(Stim.DurationMS*Nz);
% Stim.Voltage = str2num(get(handles.edit11,'string'));
% Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
% Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
% Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
% Stim.TargetRadius = 15;
% Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
% [Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
%  if get(handles.checkbox8,'Value')==1
%         handles.Setup.PointCloud.GalvoOffsetVoltage=[0,-1];
%         Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;   
%  else
%         handles.Setup.PointCloud.GalvoOffsetVoltage=[0,0];
%         Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;   
%  end
% [ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
% Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
% if get(handles.checkbox2,'Value') == 0
%     Stim.Output(:,5)=0;
% %     select = Stim.UT<2*handles.Setup.Scorepikes.sealtestduration.* Stim.UT>handles.Setup.Scorepikes.sealtestduration;
% %     Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
%     select = Stim.UT<handles.Setup.Scorepikes.sealtestduration;    
%     tempOutput=Stim.Output;
%     tempOutput(1:end,:)=0;
%     tempOutput(select,5)=1;
%     queueOutputData(handles.Setup.Daq,tempOutput);
%     startForeground(handles.Setup.Daq);
%     pause(0.1);
% else
%     Stim.Output(:,5) =0;
% end
% Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
% Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
% [Stim.LN,Stim.LX] = size(Stim.Output);
% Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
% Stim.Array = Stim.UUT-Stim.UUT;
% for i = 1:Stim.Npeaks
%     Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+i/Stim.FreqHZ);
% end
% Stim.Baseline = Stim.Output(:,1);
% if get(handles.checkbox8,'Value') == 1
%     Stim.Output(:,1)=0;
%     Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
%     temp=Stim.Output(:,2);
%     temp(temp~=0)=1;
%     Stim.Output(:,6)=temp;
% else   
%     Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
%     Stim.Output(:,2)=0;
%     temp=Stim.Output(:,1);
%     temp(temp~=0)=1;
%     Stim.Output(:,6)=temp;
% end
% Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
% Cloud.divider = floor(str2num(get(handles.edit42,'string')));
% Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
% [PXX,PYY,PZZ] = ndgrid(UX,UY,UZ);
% DataSave.PXX = PXX(:);
% DataSave.PYY = PYY(:);
% DataSave.PZZ = PZZ(:);
% 
% axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
% if ~isempty(handles.Setup.DMD.sequenceid)
%     handles.Setup=function_StopDMDSequence(handles.Setup);
% end
% for j = 1:floor(str2num(get(handles.edit9,'string'))) %do as many repetitions as needed
%     if j==1
%         for i=1:Nx*Ny
%             tic
%             if status == 0; disp('Procedure interrupted'); break; end;
%             load(['Masks11x11x21_15pixelAP\sequence' num2str(i) '.mat']);
%             handles.Setup= function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames_final))*255);
%             toc;
%         end
%     end
%     
%     handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
%     for i=1:Nx*Ny
%         if status == 0; disp('Procedure interrupted'); break; end;
%         Currentsequenceid=handles.Setup.DMD.sequenceid(i);
%         handles.Setup = function_StartProj_DMD(handles.Setup, Currentsequenceid,0);
%         queueOutputData(handles.Setup.Daq,Stim.Output);
%         Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pinout
%         if get(handles.checkbox1,'Value') == 1
%             filename = cell2mat([handles.Setup.SavingPath, handles.Setup.DataFolder, 'XYZ15AP_Seq', num2str(i), '_Rep', num2str(j), '_task', get(handles.edit46,'string'), '.mat']);
%             save(filename,'Data');
%             disp(['Data point ' num2str(i) ', ' num2str(j) ' Saved']);
%         end
%         handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
%     end
%     disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
% end
% try
% ToSave.type = 'Digital PPSF Full XZ plane';
% ToSave.Stim = Stim;
% ToSave.status = status;
% ToSave.Data = DataSave;
% catch; 
% end;
% if get(handles.checkbox1,'Value') == 1
%     filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
%     save(filename,'ToSave');
%     disp('Data Saved')
% end
% disp('Plot Result...');
% load('MaskCoordinates/Ind.mat');
% seqL=21;%each sequence contains masks for seqL points
% Score4D=zeros(Nx,Ny,Nz,str2num(get(handles.edit9,'string')));
% for j=1:str2num(get(handles.edit9,'string'))
%     if status == 0; disp('Procedure interrupted'); break; end;
%     for i=1:Nx*Ny
%         if status == 0; disp('Procedure interrupted'); break; end;
%         SavedData=load(cell2mat([handles.Setup.SavingPath, handles.Setup.DataFolder, 'XYZ15AP_Seq', num2str(i), '_Rep',num2str(j), '_task',  get(handles.edit46,'string'), '.mat']));%file name: Data
%         Data_temp=SavedData.Data.*Stim.Array';
%         Stim.baseline=mean(SavedData.Data(1:50));
%         Data_temp(Data_temp==0)=[];
%         DataMatrix=reshape(Data_temp,seqL,numel(Data_temp)/seqL);
%         for k=1:seqL
%         if status == 0; disp('Procedure interrupted'); break; end;
%         [ Score4D(Ind(k+(i-1)*seqL,1),Ind(k+(i-1)*seqL,2),Ind(k+(i-1)*seqL,3),j), ~] = function_scorespikes( handles.Setup,Stim,DataMatrix(k,:));
%         end
%     end
%     axes(handles.axes2); imagesc(UY, UX, mean(squeeze(Score4D(:,:,(Nz+1)/2,1:j)),3)); xlabel('X pixels'); ylabel('Y pixels'); title(['XY PPSF mean 1:' num2str(j)]);
%     axes(handles.axes5); imagesc(UZ, UX, mean(squeeze(Score4D(:,(Ny+1)/2,:,1:j)),3)); xlabel('X pixels'); ylabel('Z pixels'); title(['XZ PPSF mean 1:' num2str(j)]);
% end
% if get(handles.checkbox1,'Value') == 1
%     filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_Score4D_', int2str(SaveID), '_.mat'];
%     disp('Score4D Saved')
% end
% 
% end


% --- 3D digital PPSF load patterns
function pushbutton35_Callback(hObject, eventdata, handles)
global DAQstate;
global status;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
status = 1;
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);

% ProjMode = function_CheckProjMode_DMD(handles.Setup);
% if ProjMode==2301
%     [handles.Setup]=function_StopProj_DMD(handles.Setup);
%     [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
% end

Preload.X=load(['MaskCoordinates\UX.mat']);
Preload.Y=load(['MaskCoordinates\UY.mat']);
Preload.Z=load(['MaskCoordinates\UZ.mat']);
seqL=21;
if rem(numel(Preload.X.UX),seqL)~=0
    seqL=1;
end
NumofSeq=numel(Preload.X.UX)/seqL;

Cloud.Preload.sequenceid=zeros(NumofSeq,1);
for k=1:NumofSeq
    if status == 0
        disp('Procedure interrupted'); break; 
    elseif status ==2
        f = figure;set(f,'Position',[800, 800, 400,200]);
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                      'Callback','uiresume(gcbf)');
        uiwait(gcf); 
        status = 1;
        close(f);
    end
    DMDFrames_final = gpuArray(false(handles.Setup.DMD.LX,handles.Setup.DMD.LY,...
    handles.Setup.PointCloud.divider*seqL));
        for i=1:seqL
            if status == 0
                disp('Procedure interrupted'); break; 
            elseif status ==2
                f = figure;set(f,'Position',[800, 800, 400,200]);
                h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                              'Callback','uiresume(gcbf)');
                uiwait(gcf); 
                status = 1;
                close(f);
            end
                DMDFrames = function_makespots(handles.Setup,Preload.X.UX(i+(k-1)*seqL),Preload.Y.UY(i+(k-1)*seqL),Preload.Z.UZ(i+(k-1)*seqL),str2num(get(handles.edit19,'string')));   
                DMDFrames_final(:,:,((i-1)*handles.Setup.PointCloud.divider+1):(i*handles.Setup.PointCloud.divider))=DMDFrames;
        end
            [handles.Setup,Cloud.Preload.sequenceid(k)]= function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames_final))*255);
end
    disp('Random patterns are loaded!');
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
fclose(handles.Setup.SutterStage);
function_StopDMDSequence(handles.Setup);
function_Stop_DMD(handles.Setup);
delete(hObject);
end


% --- Executes random 3D PPSF.
function pushbutton36_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'DMD Digital PPSF Random';
outputSingleScan(handles.Setup.Daq,DAQstate);
% set(handles.slider3,'Value',0);
% set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);

ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

Preload.X=load(['MaskCoordinates\UX.mat']);
Preload.Y=load(['MaskCoordinates\UY.mat']);
Preload.Z=load(['MaskCoordinates\UZ.mat']);

UX = unique(Preload.X.UX);
UY = unique(Preload.Y.UY);
UZ = unique(Preload.Z.UZ);
% UX = linspace(-50,50,11);
% UY = linspace(-50,50,11);
% UZ = linspace(-100,100,21);

DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;
Nx=numel(UX);Ny=numel(UY);Nz=numel(UZ);

seqL=21;
% if rem(numel(Preload.X.UX),seqL)~=0
%     seqL=1;
% end
NumofSeq=numel(Preload.X.UX)/seqL;

Stim.Npeaks = seqL;% instead of repeat one holograph, change holograph equals to the sequence length
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
Stim.CropMask=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.CropMask = Stim.CropMask+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+i/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.nonzerovalues = find(Stim.Output(:,2));
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
    Stim.nonzerovalues = find(Stim.Output(:,1));
end
Stim.Output(:,5)=0;
select = Stim.UUT<handles.Setup.Scorepikes.sealtestduration;
Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
[PXX,PYY,PZZ] = ndgrid(UX,UY,UZ);
DataSave.PXX = PXX(:);
DataSave.PYY = PYY(:);
DataSave.PZZ = PZZ(:);

DataSave.XYZ = {}; 
load('MaskCoordinates/Ind.mat');
Score4D=zeros(Nx,Ny,Nz,str2num(get(handles.edit9,'string')));
handles.Setup.DMD.SequenceControl.RepeatModeValue=1;

axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
for j = 1:floor(str2num(get(handles.edit9,'string'))) %do as many repetitions as needed   
    handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
    for i=1:NumofSeq
        if status == 0
            disp('Procedure interrupted'); break; 
        elseif status ==2
            f = figure;set(f,'Position',[800, 800, 400,200]);
            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                          'Callback','uiresume(gcbf)');
            uiwait(gcf); 
            status = 1;
            close(f);
        end
        Currentsequenceid=Cloud.Preload.sequenceid(i);
        handles.Setup = function_StartProj_DMD(handles.Setup, Currentsequenceid);
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
        DataSave.XYZ{i,j}=Data;
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
        axes(handles.axes5);plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage');title([num2str(seqL) ' spot measurements']);
        disp(['Finish #' num2str(i) '/' num2str(NumofSeq) ' of repeat #' num2str(j)]);
    end
    disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
end
    for i=1:NumofSeq
        if status == 0
            disp('Procedure interrupted'); break; 
        elseif status ==2
            f = figure;set(f,'Position',[800, 800, 400,200]);
            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                          'Callback','uiresume(gcbf)');
            uiwait(gcf); 
            status = 1;
            close(f);
        end
        Data_temp=DataSave.XYZ{i,j}.*Stim.CropMask';
        Stim.baseline=mean(DataSave.XYZ{i,j}(1:50));
        Data_temp(Data_temp==0)=[];
        if rem(numel(Data_temp),seqL)~=0
           Data_temp=cat(1,Data_temp,0);
        end
        DataMatrix=reshape(Data_temp,numel(Data_temp)/seqL,seqL);
        DataMatrix=DataMatrix';
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        for k=1:seqL
        [ Score4D(Ind(k+(i-1)*seqL,1),Ind(k+(i-1)*seqL,2),Ind(k+(i-1)*seqL,3),j), ~] = function_scorespikes( handles.Setup,Stim,DataMatrix(k,:));
        end
    end
    axes(handles.axes2); imagesc(UY, UX, mean(squeeze(Score4D(:,:,(Nz+1)/2,1:j)),3)); xlabel('X pixels'); ylabel('Y pixels'); title(['XY PPSF mean 1:' num2str(j)]);
    axes(handles.axes3); imagesc(UZ, UX, mean(squeeze(Score4D(:,(Ny+1)/2,:,1:j)),3)); xlabel('Z pixels'); ylabel('X pixels'); title(['XZ PPSF mean 1:' num2str(j)]);

try
ToSave.type = 'Digital PPSF Full XYZ Volume random';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
ToSave.Score = Score4D;
figure(100);set(gcf, 'Position',  [100, 100, 500, 400]);imshow3D( mean(Score4D,4) );
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
catch; 
end;

% ToSave.Score = Score4D;
% figure(200);imshow3D(mean(Score4D,4));
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);

if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end

end


% --- Digital 3D PPSF XYZ sequential.
function pushbutton37_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'DMD Digital 3D PPSF Sequence';
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;
Nx=floor(str2num(get(handles.edit23,'string')));
Ny=floor(str2num(get(handles.edit24,'string')));
Nz=floor(str2num(get(handles.edit25,'string')));

Stim.Npeaks = Nz; %force speak to be 1
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
Stim.CropMask=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.CropMask = Stim.CropMask+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+i/Stim.FreqHZ-Stim.DurationMS/10000);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.nonzerovalues = find(Stim.Output(:,2));
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
    Stim.nonzerovalues = find(Stim.Output(:,1));
end
Stim.Output(:,5)=0;
select = Stim.UUT<handles.Setup.Scorepikes.sealtestduration;
Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
[PXX,PYY,PZZ] = ndgrid(UX,UY,UZ);
DataSave.PXX = PXX(:);
DataSave.PYY = PYY(:);
DataSave.PZZ = PZZ(:);
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

DataSave.XYZ = {};
[xind,yind,zind]=ind2sub([Nx,Ny,Nz],(1:Nx*Ny*Nz)');
Ind=[xind,yind,zind];
Score4D=zeros(Nx,Ny,Nz,str2num(get(handles.edit9,'string')));

axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
for j = 1:floor(str2num(get(handles.edit9,'string'))) %do as many repetitions as needed
    if j==1
        Cloud.sequenceid3D=zeros(Nx*Ny*Nz,1);
            for xx=1:Nx
                if status == 0
                    disp('Procedure interrupted'); break; 
                elseif status ==2
                    f = figure;set(f,'Position',[800, 800, 400,200]);
                    h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                                  'Callback','uiresume(gcbf)');
                    uiwait(gcf); 
                    status = 1;
                    close(f);
                end
                for yy=1:Ny
                    if status == 0
                        disp('Procedure interrupted'); break; 
                    elseif status ==2
                        f = figure;set(f,'Position',[800, 800, 400,200]);
                        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                                      'Callback','uiresume(gcbf)');
                        uiwait(gcf); 
                        status = 1;
                        close(f);
                    end
                        DMDFrames_final = gpuArray(false(handles.Setup.DMD.LX,handles.Setup.DMD.LY,...
                        handles.Setup.PointCloud.divider*Nz));
                    for i=1:Nz
                        if status == 0
                            disp('Procedure interrupted'); break; 
                        elseif status ==2
                            f = figure;set(f,'Position',[800, 800, 400,200]);
                            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                                          'Callback','uiresume(gcbf)');
                            uiwait(gcf); 
                            status = 1;
                            close(f);
                        end
                        DMDFrames = function_makespots(handles.Setup,str2num(get(handles.edit16,'string'))+UX(xx),str2num(get(handles.edit17,'string'))+UY(yy),str2num(get(handles.edit18,'string'))+UZ(i),str2num(get(handles.edit19,'string')));   
                        DMDFrames_final(:,:,((i-1)*handles.Setup.PointCloud.divider+1):(i*handles.Setup.PointCloud.divider))=DMDFrames;
                    end
                    [handles.Setup,handles.sequenceid3D((xx-1)*Nx+yy)] = function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames_final))*255);
                end
            end
    end
    for i=1:Nx*Ny
        if status == 0
            disp('Procedure interrupted'); break; 
        elseif status ==2
            f = figure;set(f,'Position',[800, 800, 400,200]);
            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                          'Callback','uiresume(gcbf)');
            uiwait(gcf); 
            status = 1;
            close(f);
        end
        handles.Setup = function_StartProj_DMD(handles.Setup, handles.sequenceid3D(i));
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
%         Data1=Data.*Stim.Array';
%         Data1(Data1==0)=[];
%         DataMatrix=reshape(Data1,Nz,numel(Data1)/Nz);
%         axes(handles.axes5);imagesc(DataMatrix);xlabel('Time [pixel]'); ylabel('Z pixels');
        DataSave.XYZ{i,j}=Data;
        axes(handles.axes5);plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage');title([num2str(i) ' coordinate']);
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
    end
    
    for i=1:Nx*Ny
        Data_temp=DataSave.XYZ{i,j}.*Stim.CropMask';
        Stim.baseline=mean(DataSave.XYZ{i,j}(1:100));
        Data_temp(Data_temp==0)=[];
        if rem(numel(Data_temp),Nz)~=0
           Data_temp=cat(1,Data_temp,0);
        end
        DataMatrix=reshape(Data_temp,numel(Data_temp)/Nz,Nz);
        DataMatrix=DataMatrix';
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        for k=1:Nz
            [ Score4D(Ind(k+(i-1)*Nz,1),Ind(k+(i-1)*Nz,2),Ind(k+(i-1)*Nz,3),j), ~] = function_scorespikes( handles.Setup,Stim,DataMatrix(k,:));
        end
    end
    axes(handles.axes2); imagesc(UY, UX, mean(squeeze(Score4D(:,:,(Nz+1)/2,1:j)),3)); xlabel('X pixels'); ylabel('Y pixels'); title(['XY PPSF mean 1:' num2str(j)]);
    axes(handles.axes3); imagesc(UZ, UX, mean(squeeze(Score4D(:,(Ny+1)/2,:,1:j)),3)); xlabel('Z pixels'); ylabel('X pixels'); title(['XZ PPSF mean 1:' num2str(j)]);
    disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
end
try
ToSave.type = 'Digital PPSF Full XYZ Volume sequential';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
ToSave.Score = Score4D;
figure(300);set(gcf, 'Position',  [100, 100, 500, 400]);imshow3D( mean(Score4D,4) );
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
catch; 
end;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end
end



function edit47_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit48_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit49_Callback(hObject, eventdata, handles)
end

function edit49_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Full XY 2D mapping at z=0
function pushbutton42_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'DMD Digital 2D XY PPSF Sequence';
ToSave.Score=[];
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))))/16;
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))))/16;
end
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;
Nx=floor(str2num(get(handles.edit23,'string')));
Ny=floor(str2num(get(handles.edit24,'string')));
Nz=floor(str2num(get(handles.edit25,'string')));
Stim.repeatNum=floor(str2num(get(handles.edit9,'string')));
[PXX,PYY,PZZ] = ndgrid(UX,UY,UZ);
DataSave.PXX = PXX(:);
DataSave.PYY = PYY(:);
DataSave.PZZ = PZZ(:);
DataSave.XY={};
% [xind,yind]=ind2sub([Nx,Ny],(1:Nx*Ny)');%sequential
% [xind,yind]=ind2sub([Nx,Ny],(randperm(Nx*Ny))');%random order
% Ind=[xind,yind];
Stim.NumofSpot=str2num(get(handles.edit57,'string'));
Stim.NumofMask=floor(Nx*Ny/Stim.NumofSpot);
% if Nx>=20
%     Ind = function_RandomIndexPoisson(Nx,Ny);
%     xind=Ind(:,1);yind=Ind(:,2);
    [ytemp1,xtemp1]=meshgrid(1:Nx,1:Ny);
    xind=xtemp1(:);
    yind=ytemp1(:);
    Ind=[xind,yind];
% else
%     [xind,yind]=ind2sub([Nx,Ny],(randperm(Nx*Ny))');%random order
%     Ind=[xind,yind];
% end

%build individual masks and store it as a sparse matrix
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
temp=((Setup.XX).^2+(Setup.YY).^2)<str2double(get(handles.edit19,'string'))^2;
DMDFrames = spalloc(handles.Setup.DMD.LX*handles.Setup.DMD.LY,Nx*Ny,nnz(temp)*Nx*Ny);
disp(['Start calculating...will take about ' num2str(0.016*Nx*Ny) 's']);
for i=1:Nx*Ny
    CenterX=str2double(get(handles.edit16,'string'))+UX(xind(i));
    CenterY=str2double(get(handles.edit17,'string'))+UY(yind(i));
    temp=((Setup.XX-CenterX).^2+(Setup.YY-CenterY).^2)<str2double(get(handles.edit19,'string'))^2;
    DMDFrames(:,i)=sparse(temp(:));
end
Score3D=zeros(Nx,Ny,Nz);

if Stim.NumofMask>450 % limit of matlab ram
    temp1=1:ceil(Stim.NumofMask/2);
    temp2=temp1(rem(Stim.NumofMask,temp1)==0);
    temp3=400-temp2;
    temp3(temp3<0)=[];
    seqL=400-min(temp3);
    NumofSeq=Stim.NumofMask/seqL;
    Stim.Npeaks=seqL;
else
    NumofSeq=1;
    seqL=Stim.NumofMask;
    Stim.Npeaks = Stim.NumofMask;
end
        
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
if numel(Stim.Voltageramp)>8
    Stim.Voltageramp=Stim.Voltage;
end
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end

[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
% figure();plot(Stim.Output(:,1));title('check whether it is 10 rising edges');
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
Stim.CropMask=Stim.Array;
% PowerMatrix=sqrt(PXX.^2+PYY.^2)/sqrt(max(UX)^2+max(UY)^2)*(str2num(get(handles.edit50,'string'))-Stim.Voltage)+Stim.Voltage;
Stim.PowerMask=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.CropMask = Stim.CropMask+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+i/Stim.FreqHZ-Stim.DurationMS/10000);
%     Stim.PowerMask=Stim.PowerMask+PowerMatrix(i)*double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
% if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
%     Stim.Output(:,1)=0;
%     Stim.Output(:,2) = Stim.PowerMask'.*Stim.Array';
%     Stim.nonzerovalues = find(Stim.Output(:,2));
% else   
%     Stim.Output(:,1) = Stim.PowerMask'.*Stim.Array';
%     Stim.Output(:,2)=0;
%     Stim.nonzerovalues = find(Stim.Output(:,1));
% end
Stim.Output(:,5)=0;
select = Stim.UUT<handles.Setup.Scorepikes.sealtestduration;
Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
Stim.Output(:,6)=Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
Cloud.XYZmap=zeros(Nx,Ny,Nz,str2num(get(handles.edit9,'string')));
Cloud.ScoreMap=zeros(Nx,Ny,Nz);

ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

Stim.CurrentXYZ = getPosition(handles.Setup.SutterStage);
disp(['Current stage Z position:' num2str(Stim.CurrentXYZ(3))]);
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes5); hold off; cla;
axes(handles.axes4); hold off; cla;
Stim.IndSeq=zeros(Nx*Ny,floor(str2num(get(handles.edit9,'string'))));

if get(handles.checkbox15,'Value') == 1
   IndPower = zeros(Nx*Ny, numel(Stim.Voltageramp));
   for i = 1:Nx*Ny
       rng('shuffle');
       IndPower(i,:)=randperm(numel(Stim.Voltageramp));
   end
   Stim.IndPower = IndPower(:);
   Stim.PowerOutput = zeros(numel(Stim.UUT),NumofSeq*numel(Stim.Voltageramp));
   ScoreXYV = zeros(Nx, Ny, numel(Stim.Voltageramp));
   for j = 1:NumofSeq*numel(Stim.Voltageramp)
       temp1 = zeros(size(Stim.Array));
       for i = 1:Stim.Npeaks
         temp1 = temp1+Stim.Voltageramp(Stim.IndPower(i+(j-1)*Stim.Npeaks))*double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
       end
       Stim.PowerOutput(:,j) = temp1';
   end
    
    % generate DMD masks
   for k=1:NumofSeq
        if status == 0
            disp('Procedure interrupted'); break; 
        elseif status ==2
            f = figure;set(f,'Position',[800, 800, 400,200]);
            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                          'Callback','uiresume(gcbf)');
            uiwait(gcf); 
            status = 1;
            close(f);
        end
        temp2 = function_RandomIndexPoisson(Nx,Ny);
        Stim.IndSeq=sub2ind([Nx,Ny],temp2(:,1),temp2(:,2));
        DMDFrames_final = gpuArray(false(handles.Setup.DMD.LX,handles.Setup.DMD.LY,seqL));
        for i=1:seqL
            if status == 0
                disp('Procedure interrupted'); break; 
            elseif status ==2
                f = figure;set(f,'Position',[800, 800, 400,200]);
                h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                              'Callback','uiresume(gcbf)');
                uiwait(gcf); 
                status = 1;
                close(f);
            end
            tempInd=Stim.IndSeq((1+(k-1)*seqL+(i-1)*Stim.NumofSpot):((k-1)*seqL+i*Stim.NumofSpot));
            DMDFrames_final(:,:,i)=full(reshape(sum(DMDFrames(:,tempInd),2),handles.Setup.DMD.LX,handles.Setup.DMD.LY));
        end
        [handles.Setup,handles.sequenceid2DXY(k)] = function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames_final))*255);
        disp(['--loading mask: ' num2str(k) '/' num2str(NumofSeq) ' finished']);
   end
for j = 1:Stim.repeatNum %do as many repetitions as needed
   for nz = 1:Nz
        if status == 0; disp('Procedure interrupted'); break; end;
        moveTime = moveTo(handles.Setup.SutterStage,[Stim.CurrentXYZ(1);Stim.CurrentXYZ(2);Stim.CurrentXYZ(3)+UZ(nz)]);
        pause(0.1);
        for v=1:numel(Stim.Voltageramp)
           for i=1:NumofSeq
                if status == 0
                    disp('Procedure interrupted'); break; 
                elseif status ==2
                    f = figure;set(f,'Position',[800, 800, 400,200]);
                    h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                                  'Callback','uiresume(gcbf)');
                    uiwait(gcf); 
                    status = 1;
                    close(f);
                end
                if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
                   Stim.Output(:,1)=0;
                   Stim.Output(:,2) = Stim.PowerOutput(:,i+(v-1)*NumofSeq);
                   Stim.nonzerovalues = find(Stim.Output(:,2));
                else   
                   Stim.Output(:,1) = Stim.PowerOutput(:,i+(v-1)*NumofSeq);
                   Stim.Output(:,2)=0;
                   Stim.nonzerovalues = find(Stim.Output(:,1));
                end
                handles.Setup = function_StartProj_DMD(handles.Setup, handles.sequenceid2DXY(i));
                queueOutputData(handles.Setup.Daq,Stim.Output);
                Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
                DataSave.XY{i+(v-1)*NumofSeq,j,nz}=Data;
                axes(handles.axes5);plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage');title([num2str(seqL) ' coordinate']);
                handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
                disp(['sequence #' num2str(i) '/' num2str(NumofSeq) ' finished!']);pause(0.5);
           end
           disp(['random permutation voltage #' num2str(v) '/' num2str(numel(Stim.Voltageramp)) ' finished!']);pause(0.5);
        end
   %plot the result
   for v=1:numel(Stim.Voltageramp)
       for i=1:NumofSeq
            if status == 0
                disp('Procedure interrupted'); break; 
            elseif status ==2
                f = figure;set(f,'Position',[800, 800, 400,200]);
                h = uicontrol('Position',[20 20 200 40],'String','Continue','Callback','uiresume(gcbf)');
                uiwait(gcf); 
                status = 1;
                close(f);
            end
            Data_temp=DataSave.XY{i+(v-1)*NumofSeq,j,nz};
            DataSave.sealtest_result(i+(v-1)*NumofSeq+(nz-1)*NumofSeq*numel(Stim.Voltageramp)+(j-1)*NumofSeq*numel(Stim.Voltageramp)*Nz) = max(Data_temp(1:2600));
            Data_temp=medfilt1(Data_temp,10);
            Data_temp=Data_temp.*Stim.CropMask';
            Data_temp(Data_temp==0)=[];
            DataMatrix=reshape(Data_temp,numel(Data_temp)/seqL,seqL);
            DataMatrix=DataMatrix';
            if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;

            for k=1:seqL
                if get(handles.checkbox14,'Value') == 0
                    tempInd=Stim.IndSeq((1+(i-1)*seqL+(k-1)*Stim.NumofSpot):((i-1)*seqL+k*Stim.NumofSpot));
                else
                    tempInd=Stim.IndSeq((1+(i-1)*seqL+(k-1)*Stim.NumofSpot):((i-1)*seqL+k*Stim.NumofSpot),j);
                end
                for n=1:Stim.NumofSpot
                  [ ScoreXYV(xind(tempInd(n)),yind(tempInd(n)),Stim.IndPower(k+(i+(v-1)*NumofSeq-1)*seqL)), ~] = function_scorespikes( handles.Setup,Stim,DataMatrix(k,:));
                end
            end
       end
   end
    figure(100);set(gcf, 'Position',  [100, 100, 300*Nz, numel(Stim.Voltageramp)*300]);
    for v = 1: numel(Stim.Voltageramp)
        subplot(numel(Stim.Voltageramp),Nz,(nz+Nz*(v-1)));
        imagesc(UX, UY, log(ScoreXYV(:,:,v)));colorbar;axis image;
        xlabel('X pixels'); ylabel('Y pixels'); title(['Z' num2str(UZ(nz)*16) '\mum, #' num2str(j), ' V=' num2str(Stim.Voltageramp(v),2)]);
        Cloud.XYZmap(:,:,nz,j,v) = ScoreXYV(:,:,v);
    end
    ToSave.Score{nz, j}=ScoreXYV;
    S1=max(squeeze(mean(Cloud.XYZmap(:,:,nz,:,:),4)),[],3);
    axes(handles.axes2); imagesc(UX, UY, log(S1)); xlabel('X pixels'); ylabel('Y pixels'); title(['Average MaxPowerProjection 1:' num2str(j)]);colorbar;
    axes(handles.axes4); plot(DataSave.sealtest_result);ylabel('nA');title('Cell Vital Signs Monitor');
    disp(['Z plane #' num2str(nz) '/' num2str(Nz) ' finished!']);
   end
    S=max(squeeze(mean(Cloud.XYZmap(:,:,:,:,end),4)),[],3);
    axes(handles.axes3); imagesc(UX, UY, log(S)); xlabel('X pixels'); ylabel('Y pixels'); title(['Average MZP at power' num2str(Stim.Voltageramp(v))]);colorbar;
    disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
end
    
else
for v=1:numel(Stim.Voltageramp)
    if status == 0; disp('Procedure interrupted'); break; end;
    if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
       Stim.Output(:,1)=0;
       Stim.Output(:,2) = Stim.Voltageramp(v).*Stim.Array';
       Stim.nonzerovalues = find(Stim.Output(:,2));
    else   
       Stim.Output(:,1) = Stim.Voltageramp(v).*Stim.Array';
       Stim.Output(:,2)=0;
       Stim.nonzerovalues = find(Stim.Output(:,1));
    end
for j = 1:Stim.repeatNum %do as many repetitions as needed
    if status == 0; disp('Procedure interrupted'); break; end;
    for nz = 1:Nz
        if status == 0; disp('Procedure interrupted'); break; end;
        moveTime = moveTo(handles.Setup.SutterStage,[Stim.CurrentXYZ(1);Stim.CurrentXYZ(2);Stim.CurrentXYZ(3)+UZ(nz)]);
        pause(0.1);
    if get(handles.checkbox14,'Value') == 0    
        if j==1 && nz==1 && v==1
            for k=1:NumofSeq
            if status == 0
                disp('Procedure interrupted'); break; 
            elseif status ==2
                f = figure;set(f,'Position',[800, 800, 400,200]);
                h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                              'Callback','uiresume(gcbf)');
                uiwait(gcf); 
                status = 1;
                close(f);
            end
            temp2 = function_RandomIndexPoisson(Nx,Ny);
            Stim.IndSeq=sub2ind([Nx,Ny],temp2(:,1),temp2(:,2));
            DMDFrames_final = gpuArray(false(handles.Setup.DMD.LX,handles.Setup.DMD.LY,seqL));
                for i=1:seqL
                    if status == 0
                        disp('Procedure interrupted'); break; 
                    elseif status ==2
                        f = figure;set(f,'Position',[800, 800, 400,200]);
                        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                                      'Callback','uiresume(gcbf)');
                        uiwait(gcf); 
                        status = 1;
                        close(f);
                    end
                    tempInd=Stim.IndSeq((1+(k-1)*seqL+(i-1)*Stim.NumofSpot):((k-1)*seqL+i*Stim.NumofSpot));
                    DMDFrames_final(:,:,i)=full(reshape(sum(DMDFrames(:,tempInd),2),handles.Setup.DMD.LX,handles.Setup.DMD.LY));
                end
                [handles.Setup,handles.sequenceid2DXY(k)] = function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames_final))*255);
                disp(['--loading mask: ' num2str(k) '/' num2str(NumofSeq) ' finished']);
            end
        end
    else
        for k=1:NumofSeq
            if status == 0
                disp('Procedure interrupted'); break; 
            elseif status ==2
                f = figure;set(f,'Position',[800, 800, 400,200]);
                h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                              'Callback','uiresume(gcbf)');
                uiwait(gcf); 
                status = 1;
                close(f);
            end
            temp2 = function_RandomIndexPoisson(Nx,Ny);
            Stim.IndSeq(:,j)=sub2ind([Nx,Ny],temp2(:,1),temp2(:,2));
            DMDFrames_final = false(handles.Setup.DMD.LX,handles.Setup.DMD.LY,seqL);
            for i=1:seqL
                tempInd=Stim.IndSeq((1+(k-1)*seqL+(i-1)*Stim.NumofSpot):((k-1)*seqL+i*Stim.NumofSpot),j);
                DMDFrames_final(:,:,i)=full(reshape(sum(DMDFrames(:,tempInd),2),handles.Setup.DMD.LX,handles.Setup.DMD.LY));
            end
            [handles.Setup,handles.sequenceid2DXY(k)] = function_StoreImages_DMD(handles.Setup, uint8(DMDFrames_final)*255);
            disp(['--loading mask: ' num2str(k) '/' num2str(NumofSeq) ' finished']);
        end
    end
    
    for i=1:NumofSeq
        if status == 0
            disp('Procedure interrupted'); break; 
        elseif status ==2
            f = figure;set(f,'Position',[800, 800, 400,200]);
            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                          'Callback','uiresume(gcbf)');
            uiwait(gcf); 
            status = 1;
            close(f);
        end
        handles.Setup = function_StartProj_DMD(handles.Setup, handles.sequenceid2DXY(i));
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
        DataSave.XY{i,j,nz,v}=Data;
        axes(handles.axes5);plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage');title([num2str(seqL) ' coordinate']);
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
        disp(['sequence #' num2str(i) '/' num2str(NumofSeq) ' finished!']);pause(1);
    end
    
    for i=1:NumofSeq
        if status == 0
            disp('Procedure interrupted'); break; 
        elseif status ==2
            f = figure;set(f,'Position',[800, 800, 400,200]);
            h = uicontrol('Position',[20 20 200 40],'String','Continue','Callback','uiresume(gcbf)');
            uiwait(gcf); 
            status = 1;
            close(f);
        end
        Data_temp=DataSave.XY{i,j,nz,v};
        DataSave.sealtest_result(i+(nz-1)*NumofSeq+(j-1)*NumofSeq*Nz+(v-1)*NumofSeq*Nz*Stim.repeatNum) = max(Data_temp(1:2600));
        Data_temp=medfilt1(Data_temp,10);
        Data_temp=Data_temp.*Stim.CropMask';
        Data_temp(Data_temp==0)=[];
        DataMatrix=reshape(Data_temp,numel(Data_temp)/seqL,seqL);
        DataMatrix=DataMatrix';
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        for k=1:seqL
            if get(handles.checkbox14,'Value') == 0
                tempInd=Stim.IndSeq((1+(i-1)*seqL+(k-1)*Stim.NumofSpot):((i-1)*seqL+k*Stim.NumofSpot));
            else
                tempInd=Stim.IndSeq((1+(i-1)*seqL+(k-1)*Stim.NumofSpot):((i-1)*seqL+k*Stim.NumofSpot),j);
            end
            for n=1:Stim.NumofSpot
              [ Score3D(xind(tempInd(n)),yind(tempInd(n)),nz), ~] = function_scorespikes( handles.Setup,Stim,DataMatrix(k,:));
            end
        end
    end
    figure(100);set(gcf, 'Position',  [100, 100, 300*Nz, numel(Stim.Voltageramp)*300]);
    subplot(numel(Stim.Voltageramp),Nz,(nz+Nz*(v-1)));
    imagesc(UX, UY, log(Score3D(:,:,nz)));colorbar;axis image;
    xlabel('X pixels'); ylabel('Y pixels'); title(['Z' num2str(UZ(nz)*16,2) '\mum, #' num2str(j), ' V=' num2str(Stim.Voltageramp(v),2)]);
end
ToSave.Score{j,v}=Score3D;
Cloud.XYZmap(:,:,:,j,v)=Score3D;
S1=max(mean(Cloud.XYZmap(:,:,:,:,v),4),[],3);
axes(handles.axes2); imagesc(UX, UY, log(S1)); xlabel('X pixels'); ylabel('Y pixels'); title(['Average MZP 1:' num2str(j)]);colorbar;
axes(handles.axes4); plot(DataSave.sealtest_result);ylabel('nA');title('Cell Vital Signs Monitor');
disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
end
S=max(mean(Cloud.XYZmap(:,:,:,:,v),4),[],3);
axes(handles.axes3); imagesc(UX, UY, log(S)); xlabel('X pixels'); ylabel('Y pixels'); title(['Average MZP at power' num2str(Stim.Voltageramp(v))]);colorbar;
disp(['---Voltage #' num2str(v) '/' num2str(numel(Stim.Voltageramp)) ' finished!']);
end
end

moveTime = moveTo(handles.Setup.SutterStage,[Stim.CurrentXYZ(1);Stim.CurrentXYZ(2);Stim.CurrentXYZ(3)]);

ScoreMap = sum(Cloud.XYZmap(:,:,:,:,numel(Stim.Voltageramp)),4);
if isempty(find(ScoreMap, 1))
   disp('No spike is found in the mapping result!');
else
    [counts,edges] = histcounts(ScoreMap);
    thres=numel(counts)-1;
    while sum(counts(thres:end))<(0.1*Nx*Ny*Nz)
        thres=thres-1;
        if thres<1
            break;
        end
    end
    ScoreMap(ScoreMap<=edges(thres+1))=0;%remove background
    Cloud.ScoreMap=ScoreMap;
end

try
ToSave.Ind=Ind;
ToSave.type = 'Digital PPSF Full XY 2D sequential';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
catch; 
end;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end
end

% --- high resolution scan sub-regions after full xy 2D scanning
function pushbutton43_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'Fine ROI mapping';
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end
Nx=floor(str2num(get(handles.edit23,'string')));
Ny=floor(str2num(get(handles.edit24,'string')));
Nz=floor(str2num(get(handles.edit25,'string')));
UZ=round(UZ/handles.Setup.DataAnalysis.Zratio);
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;

MAP=max(mean(Cloud.XYZmap,4),[],3);
[counts, bin]=hist(MAP(:),round(numel(MAP)*0.02));
Threshold=min(bin(counts<round(numel(MAP)*0.05)));%sparse signal: signal pixels should be less than 5% of total pixels.
[xind,yind]=ind2sub([Nx,Ny],find(MAP>Threshold));%use func "find", the result is sorted
CenterX=UX(xind);CenterY=UY(yind);
Cloud.CenterX=CenterX;Cloud.CenterY=CenterY;
CenterX=cat(2,CenterX,[min(CenterX)-2*Stim.TargetRadius,max(CenterX)+2*Stim.TargetRadius]);
CenterY=cat(2,CenterY,[min(CenterY)-2*Stim.TargetRadius,max(CenterY)+2*Stim.TargetRadius]);


    Num=round(sqrt(numel(CenterX))*(str2num(get(handles.edit55,'string')))^2);
    [Xq,Yq]=meshgrid(linspace(min(CenterX),max(CenterX),Num),linspace(min(CenterY),max(CenterY),Num));
    Zq=griddata(CenterX,CenterY,ones(size(CenterX)),Xq,Yq);
    Stim.Data.IND=find(~isnan(Zq));
    Stim.N=numel(Stim.Data.IND);
    if rem(Stim.N,2)==1
        Stim.Data.IND(1)=[];
        Stim.N=numel(Stim.Data.IND);
    end
    Stim.Data.IND=Stim.Data.IND(randperm(Stim.N));
    Stim.Data.CX=Xq(Stim.Data.IND);
    Stim.Data.CY=Yq(Stim.Data.IND);
    [Xq1,Yq1]=meshgrid(1:Num,1:Num);
    Ind=[Xq1(Stim.Data.IND),Yq1(Stim.Data.IND)];
    axes(handles.axes4);scatter(Stim.Data.CX,Stim.Data.CY,'r.');
    hold on;scatter(CenterX,CenterY,'bo');axis tight;
    title(['ROI, # of pixel = ' num2str(Stim.N)]);

if Stim.N>450 % limit of matlab ram
    temp1=1:ceil(Stim.N/2);
    temp2=temp1(rem(Stim.N,temp1)==0);
    temp3=400-temp2;
    temp3(temp3<0)=[];
    seqL=400-min(temp3);
    NumofSeq=Stim.N/seqL;
    Stim.Npeaks=seqL;
else
    NumofSeq=1;
    seqL=Stim.N;
    Stim.Npeaks = Stim.N;
end
        
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end

[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
Stim.CropMask=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.CropMask = Stim.CropMask+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+i/Stim.FreqHZ-Stim.DurationMS/10000);
end
Stim.Baseline = Stim.Output(:,1);
Stim.Output(:,5)=0;
select = Stim.UUT<handles.Setup.Scorepikes.sealtestduration;
Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
Stim.Output(:,6)=Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
Cloud.XYZmapVoltage=zeros(Num,Num,Nz,str2num(get(handles.edit9,'string')),numel(Stim.Voltageramp));
Score3D=zeros(Num,Num,Nz);

ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

Stim.CurrentXYZ = getPosition(handles.Setup.SutterStage);
disp(['Current stage Z position:' num2str(Stim.CurrentXYZ(3))]);
axes(handles.axes2); hold off; cla;axes(handles.axes5); hold off; cla;

for v=1:numel(Stim.Voltageramp)
    if status == 0
        disp('Procedure interrupted'); break; 
    elseif status ==2
        f = figure;set(f,'Position',[800, 800, 400,200]);
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                                  'Callback','uiresume(gcbf)');
        uiwait(gcf); 
        status = 1;
        close(f);
    end
    if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
       Stim.Output(:,1)=0;
       Stim.Output(:,2) = Stim.Voltageramp(v).*Stim.Array';
       Stim.nonzerovalues = find(Stim.Output(:,2));
    else   
       Stim.Output(:,1) = Stim.Voltageramp(v).*Stim.Array';
       Stim.Output(:,2)=0;
       Stim.nonzerovalues = find(Stim.Output(:,1));
    end
for j = 1:floor(str2num(get(handles.edit9,'string'))) %do as many repetitions as needed
    if status == 0
        disp('Procedure interrupted'); break; 
    elseif status ==2
        f = figure;set(f,'Position',[800, 800, 400,200]);
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                      'Callback','uiresume(gcbf)');
        uiwait(gcf); 
        status = 1;
        close(f);
    end
    for nz = 1:Nz
        if status == 0; disp('Procedure interrupted'); break; end;
        moveTime = moveTo(handles.Setup.SutterStage,[Stim.CurrentXYZ(1);Stim.CurrentXYZ(2);Stim.CurrentXYZ(3)+UZ(nz)]);
        pause(0.1);
    if j==1 && nz==1 && v==1
        for k=1:NumofSeq
        DMDFrames_final = gpuArray(false(handles.Setup.DMD.LX,handles.Setup.DMD.LY,seqL));
            for i=1:seqL
                if status == 0; disp('Procedure interrupted'); break; end;
                   DMDFrames = function_makespots(handles.Setup,str2num(get(handles.edit16,'string'))+Stim.Data.CX(i+(k-1)*seqL),str2num(get(handles.edit17,'string'))+Stim.Data.CY(i+(k-1)*seqL),0,str2num(get(handles.edit19,'string')));   
                   DMDFrames_final(:,:,i)=DMDFrames(:,:,1);
            end
            [handles.Setup,handles.sequenceid2DXY(k)] = function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames_final))*255);
            disp(['--loading mask: ' num2str(k) '/' num2str(NumofSeq) ' finished']);
        end
    end
    
    for i=1:NumofSeq
        if status == 0
            disp('Procedure interrupted'); break; 
        elseif status ==2
            f = figure;set(f,'Position',[800, 800, 400,200]);
            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                          'Callback','uiresume(gcbf)');
            uiwait(gcf); 
            status = 1;
            close(f);
        end
        handles.Setup = function_StartProj_DMD(handles.Setup, handles.sequenceid2DXY(i));
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
        DataSave.XY{i,j,nz,v}=Data;
        axes(handles.axes5);plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage');title([num2str(seqL) ' coordinate']);
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
        disp(['sequence #' num2str(i) '/' num2str(NumofSeq) ' finished!']);pause(1);
    end
    
    for i=1:NumofSeq
        if status == 0
            disp('Procedure interrupted'); break; 
        elseif status ==2
            f = figure;set(f,'Position',[800, 800, 400,200]);
            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                          'Callback','uiresume(gcbf)');
            uiwait(gcf); 
            status = 1;
            close(f);
        end
        Data_temp=DataSave.XY{i,j,nz,v};
        Data_temp=medfilt1(Data_temp,10);
        Data_temp=Data_temp.*Stim.CropMask';
        Data_temp(Data_temp==0)=[];
        DataMatrix=reshape(Data_temp,numel(Data_temp)/seqL,seqL);
        DataMatrix=DataMatrix';
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        for k=1:seqL
            [ Score3D(Ind(k+(i-1)*seqL,1),Ind(k+(i-1)*seqL,2),nz), ~] = function_scorespikes( handles.Setup,Stim,DataMatrix(k,:));
        end
    end
    figure(200);set(gcf, 'Position',  [100, 100, Nz*100, numel(Stim.Voltageramp)*100]);
    subplot(numel(Stim.Voltageramp),Nz,(nz+Nz*(v-1)));
    imagesc(Xq(1,:), Yq(:,1), log(Score3D(:,:,nz)'));colorbar;axis image;
    xlabel('X pixels'); ylabel('Y pixels'); title(['Z' num2str(UZ(nz),2) '\mum, #' num2str(j), ' V=' num2str(Stim.Voltageramp(v),2)]);
    Score2D=max(Score3D,[],3);
    axes(handles.axes2); imagesc(Xq(1,:), Yq(:,1), log(Score2D')); xlabel('X pixels'); ylabel('Y pixels'); title(['Maximum Z Projection:' num2str(j)]);colorbar;
end
ToSave.Score{j,v}=Score3D;
Cloud.XYZmapVoltage(:,:,:,j,v)=Score3D;
disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
end
disp(['---Voltage #' num2str(v) '/' num2str(numel(Stim.Voltageramp)) ' finished!']);
end
moveTime = moveTo(handles.Setup.SutterStage,[Stim.CurrentXYZ(1);Stim.CurrentXYZ(2);Stim.CurrentXYZ(3)]);

try
ToSave.Xq=Xq;
ToSave.Yq=Yq;
ToSave.Ind=Ind;
ToSave.type = 'Fine ROI mapping';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
catch; 
end;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved');
end
end


% --- simultaneously stimulate multiple targets in 3D
function pushbutton44_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'Stimulate multiple targets in 3D';
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;
Nx=floor(str2num(get(handles.edit23,'string')));
Ny=floor(str2num(get(handles.edit24,'string')));
Nz=floor(str2num(get(handles.edit25,'string')));

Stim.Npeaks = floor(str2num(get(handles.edit6,'string'))); % Number of blue light pulses in test
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
if get(handles.checkbox2,'Value') == 0
    Stim.Output(:,5)=0;
%     select = Stim.UT<2*handles.Setup.Scorepikes.sealtestduration.* Stim.UT>handles.Setup.Scorepikes.sealtestduration;
%     Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
    select = Stim.UT<handles.Setup.Scorepikes.sealtestduration;    
    tempOutput=Stim.Output;
    tempOutput(1:end,:)=0;
    tempOutput(select,5)=1;
    queueOutputData(handles.Setup.Daq,tempOutput);
    startForeground(handles.Setup.Daq);
    pause(0.1);
else
    Stim.Output(:,5) =0;
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value')==1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.nonzerovalues = find(Stim.Output(:,2));
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
    Stim.nonzerovalues = find(Stim.Output(:,1));
end
Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end
DataSave.MultiTargets = {};DataSave.Score={};

% MAP=mean(Cloud.XYZmap,4);
% [counts, bin]=hist(MAP(:),round(numel(MAP)*0.02));
% Threshold=min(bin(counts<round(numel(MAP)*0.02)));%sparse signal: signal pixels should be less than 5% of total pixels.
% [xind,yind,zind]=ind2sub([Nx,Ny,Nz],find(MAP>Threshold));%use func "find", the result is sorted
% CenterX=UX(xind);CenterY=UY(yind);CenterZ=UZ(zind)/handles.Setup.DataAnalysis.Zratio;
% Stim.Targets=[CenterX',CenterY',CenterZ];
% Distance=sqrt(sum(diff(Stim.Targets).^2,2));%calculate the distance between sorted spots
% Ind_temp=find(Distance<abs(UX(1)-UX(2))*4)+1;
% Stim.Targets(Ind_temp,:)=[];

Stim.Targets = xlsread('targets.xlsx',1);
Ntarget=size(Stim.Targets,1);
handles.Setup.DMD.SequenceControl.RepeatModeValue=Stim.Npeaks;
axes(handles.axes2); hold off; cla;axes(handles.axes5); hold off; cla;
% axes(handles.axes3); scatter(Stim.Targets(:,1), Stim.Targets(:,2),'blue','filled'); xlabel('X(pixel)'); ylabel('Y(pixel)'); title('Targets');
for indcomb=1:Ntarget
    Stim.IndComb{indcomb}=nchoosek(1:Ntarget,indcomb);
    for p=1:size(Stim.IndComb{indcomb},1)
        ic=Stim.IndComb{indcomb}(p,:);
        ix=Stim.Targets(ic,1);
        iy=Stim.Targets(ic,2);
        iz=Stim.Targets(ic,3);
        for j = 1:floor(str2num(get(handles.edit9,'string'))) %do as many repetitions as needed
         if status == 0; disp('Procedure interrupted'); break; end;
            [DMDFrames] =  function_makespots_ori(handles.Setup,ix,iy,iz,str2num(get(handles.edit19,'string'))*ones(size(ix)));
            handles.Setup = function_feed_DMD(handles.Setup, uint8(gather(DMDFrames))*255);
            queueOutputData(handles.Setup.Daq,Stim.Output);
            Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
            Stim.baseline=mean(Data(1:100));
            if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else handles.Setup.Scorepikes.Method=0; end;
            [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
            axes(handles.axes2);  scatter(j, score,'red','filled'); hold on; xlabel('repeat#'); ylabel('Score'); title(['Stimulate ' num2str(indcomb) ' spots']);
            axes(handles.axes5);  plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage');
            DataSave.MultiTargets{p,indcomb,j} = Data;
            DataSave.Score{p,indcomb,j} = score;
            handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
        end
        axes(handles.axes2); hold off; cla;
    end
    disp(['Combination ' num2str(indcomb) '/' num2str(Ntarget) ' finished!']);
end
try
ToSave.type = 'Simultaneous multiple 3D targets';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
catch; end;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end
end


% --- Red stimulation.
function checkbox11_Callback(hObject, eventdata, handles)
end



function edit55_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function edit55_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
global status;
status = 2; % 0: interrupt, 2: pause and show dialog window
f = figure;set(f,'Position',[800, 800, 400,200]);
h = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
uiwait(gcf); 
status = 1;
close(f);
end


% --- Fill Coarse mapping parameters.
function pushbutton46_Callback(hObject, eventdata, handles)
global DAQstate;
global status; status = 1; % 1 if ok to continue
global SaveID; SaveID=0;
global ToSave;
global Stage;
global Cloud;
global StageStatus; StageStatus=0;
DAQstate=handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.edit6,'string',int2str(1)); % Number of pulses
set(handles.edit7,'string',int2str(40)); % Frequency in Hertz
set(handles.edit8,'string',int2str(4)); % Pulse duration in ms
set(handles.edit9,'string',int2str(5)); % Number of repetitions
set(handles.edit10,'string',int2str(0)); % Delay between repeitions in seconds
set(handles.edit11,'string',num2str(3)); %Stim laser voltage,center (default, no power sweep)
set(handles.edit29,'string',num2str(0.09)); %Pre stimulation delay in seconds
set(handles.edit30,'string',num2str(0.8)); %Pulsed lase Duty cycle between 0 and 1
set(handles.edit12,'string',num2str(2.2)); %Voltage sweep start voltage
set(handles.edit13,'string',num2str(3.8));%Voltage sweep End voltage
set(handles.edit14,'string',num2str(3));%Voltage sweep number of steps
set(handles.edit16,'string',num2str(0)); %Targets positions pixels X
set(handles.edit17,'string',num2str(0)); %Targets positions pixels Y
set(handles.edit18,'string',num2str(0)); %Targets positions pixels Z
set(handles.edit19,'string',num2str(20)); %Target radius
set(handles.edit20,'string',num2str(-600)); % PPSFRange of measurements on X axis in pixels,start
set(handles.edit21,'string',num2str(-600));% PPSF Range of measurements on Y axis in pixels,start
set(handles.edit22,'string',num2str(-60));% PPSFRange of measurements on Z axis in pixels,start
set(handles.edit47,'string',num2str(600)); % PPSFRange of measurements on X axis in pixels,end
set(handles.edit48,'string',num2str(600));% PPSF Range of measurements on Y axis in pixels,end
set(handles.edit49,'string',num2str(-60));% PPSFRange of measurements on Z axis in pixels,end
set(handles.edit23,'string',num2str(20));% PPSF number of points on X axis
set(handles.edit24,'string',num2str(20));% PPSF number of points on Y axis
set(handles.edit25,'string',num2str(1));% PPSF number of points on Z axis
set(handles.edit55,'string',num2str(2.5));% Sampling grid
set(handles.checkbox1,'Value', 1);
set(handles.checkbox2,'Value', 0);
end


% --- Executes on button press in pushbutton47.
function pushbutton47_Callback(hObject, eventdata, handles)
global DAQstate;
global status; status = 1; % 1 if ok to continue
global SaveID; SaveID=0;
global ToSave;
global Stage;
global Cloud;
global StageStatus; StageStatus=0;
DAQstate=handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.edit6,'string',int2str(1)); % Number of pulses
set(handles.edit7,'string',int2str(80)); % Frequency in Hertz
set(handles.edit8,'string',int2str(4)); % Pulse duration in ms
set(handles.edit9,'string',int2str(5)); % Number of repetitions
set(handles.edit10,'string',int2str(0)); % Delay between repeitions in seconds
set(handles.edit11,'string',num2str(4)); %Stim laser voltage,center (default)
set(handles.edit29,'string',num2str(0.09)); %Pre stimulation delay in seconds
set(handles.edit30,'string',num2str(0.8)); %Pulsed lase Duty cycle between 0 and 1
set(handles.edit14,'string',num2str(4));%Voltage sweep number of steps
set(handles.edit16,'string',num2str(0)); %Targets positions pixels X
set(handles.edit17,'string',num2str(0)); %Targets positions pixels Y
set(handles.edit18,'string',num2str(0)); %Targets positions pixels Z
set(handles.edit19,'string',num2str(20)); %Target radius
set(handles.edit22,'string',num2str(-120));% PPSFRange of measurements on Z axis in pixels,start
set(handles.edit49,'string',num2str(20));% PPSFRange of measurements on Z axis in pixels,end
set(handles.edit23,'string',num2str(40));% PPSF number of points on X axis
set(handles.edit24,'string',num2str(40));% PPSF number of points on Y axis
set(handles.edit25,'string',num2str(5));% PPSF number of points on Z axis
set(handles.edit55,'string',num2str(2.5));% Sampling grid
set(handles.checkbox1,'Value', 1);
set(handles.checkbox2,'Value', 0);
end


% % --- preload patterns for random mapping
% function pushbutton50_Callback(hObject, eventdata, handles)
% global DAQstate;
% global status;
% global ToSave;
% global Cloud;
% global SaveID; SaveID=SaveID+1;
% status = 1;
% DAQstate = [0 0 0 0 0 0 0];
% outputSingleScan(handles.Setup.Daq,DAQstate);
% set(handles.slider3,'Value',0);
% set(handles.slider2,'Value',0);
% set(handles.edit3,'string',0);
% set(handles.edit4,'string',0);
% 
% % ProjMode = function_CheckProjMode_DMD(handles.Setup);
% % if ProjMode==2301
% %     [handles.Setup]=function_StopProj_DMD(handles.Setup);
% %     [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
% % end
% 
% load(['MaskCoordinates_mapping\focicases.mat']);
% load(['MaskCoordinates_mapping\seqLcases.mat']);
% Preload.X=load(['MaskCoordinates_mapping\UX.mat']);
% Preload.Y=load(['MaskCoordinates_mapping\UY.mat']);
% Preload.Z=load(['MaskCoordinates_mapping\UZ.mat']);
% NumofSpot=str2num(get(handles.edit56,'string'));
% [~,indtemp]=min(abs(focicases-NumofSpot));
% NumofSpot=focicases(indtemp);
% seqL=seqLcases(indtemp);
% if rem(numel(Preload.X.UX),seqL)~=0
%     seqL=1;
% end
% NumofSeq=numel(Preload.X.UX)/(seqL*NumofSpot);
% 
% Cloud.Preload.sequenceid=zeros(NumofSeq,1);
% disp(['Start calculating masks. It will take about ' num2str(0.4*seqL*NumofSeq) ' seconds.']);
% 
% %compute the sparse matrix of each individual spot in 3D----------------
% %----------------------
% for k=1:NumofSeq
%     if status == 0
%         disp('Procedure interrupted'); break; 
%     elseif status ==2
%         f = figure;set(f,'Position',[800, 800, 400,200]);
%         h = uicontrol('Position',[20 20 200 40],'String','Continue',...
%                       'Callback','uiresume(gcbf)');
%         uiwait(gcf); 
%         status = 1;
%         close(f);
%     end
%     DMDFrames_final = gpuArray(false(handles.Setup.DMD.LX,handles.Setup.DMD.LY,...
%     handles.Setup.PointCloud.divider*seqL));
%         for i=1:seqL
%             if status == 0
%                 disp('Procedure interrupted'); break; 
%             elseif status ==2
%                 f = figure;set(f,'Position',[800, 800, 400,200]);
%                 h = uicontrol('Position',[20 20 200 40],'String','Continue',...
%                               'Callback','uiresume(gcbf)');
%                 uiwait(gcf); 
%                 status = 1;
%                 close(f);
%             end
%                 tempx=Preload.X.UX((1+(k-1)*seqL+(i-1)*NumofSpot):((k-1)*seqL+i*NumofSpot));
%                 tempy=Preload.Y.UY((1+(k-1)*seqL+(i-1)*NumofSpot):((k-1)*seqL+i*NumofSpot));
%                 tempz=Preload.Z.UZ((1+(k-1)*seqL+(i-1)*NumofSpot):((k-1)*seqL+i*NumofSpot));
%                 [DMDFrames] = function_makespots_ori(handles.Setup,tempx,tempy,tempz,ones(size(tempx))*str2num(get(handles.edit19,'string')));
%                 DMDFrames_final(:,:,((i-1)*handles.Setup.PointCloud.divider+1):(i*handles.Setup.PointCloud.divider))=DMDFrames;
%         end
%             [handles.Setup,Cloud.Preload.sequenceid(k)]= function_StoreImages_DMD(handles.Setup, uint8(gather(DMDFrames_final))*255);
% end
%     disp('Random patterns for 4D mapping are loaded!');
% DAQstate = [0 0 0 0 0 0 0];
% outputSingleScan(handles.Setup.Daq,DAQstate);
% end

% --- preload pattern in 3D random mapping.
function pushbutton53_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'Digital random mapping 3D';
outputSingleScan(handles.Setup.Daq,DAQstate);
% set(handles.slider3,'Value',0);
% set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
Stim.repeatNum = floor(str2num(get(handles.edit9,'string')));

ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end

Nx=floor(str2num(get(handles.edit23,'string')));
Ny=floor(str2num(get(handles.edit24,'string')));
Nz=floor(str2num(get(handles.edit25,'string')));

[yytemp,xxtemp,zztemp]=meshgrid(1:Nx,1:Ny,1:Nz);
xind=xxtemp(:);
yind=yytemp(:);
zind=zztemp(:);

Stim.NumofSpot=str2num(get(handles.edit56,'string'));
Stim.NumofMask=floor(Nx*Ny*Nz/Stim.NumofSpot*handles.Setup.PointCloud.divider);
if Stim.NumofMask>450 % limit of matlab ram
    temp1=1:ceil(Stim.NumofMask/2);
    temp2=temp1(rem(Stim.NumofMask,temp1)==0);
    temp3=400-temp2;
    temp3(temp3<0)=[];
    seqL=400-min(temp3);
    NumofSeq=Stim.NumofMask/seqL;
    Stim.Npeaks=seqL;
else
    NumofSeq=1;
    seqL=Stim.NumofMask;
    Stim.Npeaks = Stim.NumofMask;
end

Cloud.IndSeq=[];
Cloud.sequenceid3DXYZ=[];
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
temp=((Setup.XX).^2+(Setup.YY).^2)<str2double(get(handles.edit19,'string'))^2;
disp(['Start calculating...will take about ' num2str(0.016*Nx*Ny*Nz*handles.Setup.PointCloud.divider) 's']);
for t=1:handles.Setup.PointCloud.divider
    DMDFrames.(sprintf('s%d',t)) = spalloc(handles.Setup.DMD.LX*handles.Setup.DMD.LY,Nx*Ny*Nz,nnz(temp)*Nx*Ny*Nz);
    for i=1:Nx*Ny*Nz
        CenterX=UX(xind(i))+UZ(zind(i))*cos(t*2*pi/handles.Setup.PointCloud.divider+handles.Setup.PointCloud.phiDMD)*0.74;%voltage =1-->0.22; voltage=2-->0.74
        CenterY=UY(yind(i))+UZ(zind(i))*sin(t*2*pi/handles.Setup.PointCloud.divider+handles.Setup.PointCloud.phiDMD)*0.74;
        temp=((Setup.XX-CenterX).^2+(Setup.YY-CenterY).^2)<str2double(get(handles.edit19,'string'))^2;
        DMDFrames.(sprintf('s%d',t))(:,i)=sparse(temp(:));
    end
    disp([num2str(t) '/' num2str(handles.Setup.PointCloud.divider) ' finished']);
end
disp('Start loading patterns...');
for j=1:Stim.repeatNum
      Cloud.IndSeq(:,j)=randperm(Nx*Ny*Nz);
        for k=1:NumofSeq
                DMDFrames_final = false(handles.Setup.DMD.LX,handles.Setup.DMD.LY,seqL);
                for i=1:seqL/handles.Setup.PointCloud.divider
                    tempInd=Cloud.IndSeq((1+(k-1)*seqL*Stim.NumofSpot/handles.Setup.PointCloud.divider+(i-1)*Stim.NumofSpot):((k-1)*seqL*Stim.NumofSpot/handles.Setup.PointCloud.divider+i*Stim.NumofSpot),j);
                    for t=1:handles.Setup.PointCloud.divider
                        DMDFrames_final(:,:,t+(i-1)*handles.Setup.PointCloud.divider)=full(reshape(sum(DMDFrames.(sprintf('s%d',t))(:,tempInd),2),handles.Setup.DMD.LX,handles.Setup.DMD.LY));
                    end
                end
                [handles.Setup,Cloud.sequenceid3DXYZ(k+(j-1)*NumofSeq)] = function_StoreImages_DMD(handles.Setup, uint8(DMDFrames_final)*255);
                disp(['--loading mask: ' num2str(k) '/' num2str(NumofSeq) ' in repeat #' num2str(j) '/' num2str(Stim.repeatNum) ' finished']);
        end
end
disp('Calculation completed!');
end

% --- Project random patterns from pre-calculated for mapping
function pushbutton51_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate =handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'Digital random mapping 3D';
outputSingleScan(handles.Setup.Daq,DAQstate);
% set(handles.slider3,'Value',0);
% set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);

ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;
Nx=floor(str2num(get(handles.edit23,'string')));
Ny=floor(str2num(get(handles.edit24,'string')));
Nz=floor(str2num(get(handles.edit25,'string')));

Stim.NumofSpot=str2num(get(handles.edit56,'string'));
Stim.NumofMask=floor(Nx*Ny*Nz/Stim.NumofSpot*handles.Setup.PointCloud.divider);
if Stim.NumofMask>450 % limit of matlab ram
    temp1=1:ceil(Stim.NumofMask/2);
    temp2=temp1(rem(Stim.NumofMask,temp1)==0);
    temp3=400-temp2;
    temp3(temp3<0)=[];
    seqL=400-min(temp3);
    NumofSeq=Stim.NumofMask/seqL;
    Stim.Npeaks=seqL;
else
    NumofSeq=1;
    seqL=Stim.NumofMask;
    Stim.Npeaks = Stim.NumofMask;
end

[PXX,PYY,PZZ] = ndgrid(UX,UY,UZ);
DataSave.PXX = PXX(:);
DataSave.PYY = PYY(:);
DataSave.PZZ = PZZ(:);

[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));   
Stim.Npeaks = seqL/handles.Setup.PointCloud.divider;% instead of repeat one holograph, change holograph equals to the sequence length
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
Stim.repeatNum=floor(str2num(get(handles.edit9,'string')));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
Stim.CropMask=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.CropMask = Stim.CropMask+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+i/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
Stim.Output(:,5)=0;
select = Stim.UUT<handles.Setup.Scorepikes.sealtestduration;
Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));

DataSave.XYZ = {}; 
% load('MaskCoordinates_mapping/Ind.mat');
Score4D=zeros(Nx,Ny,Nz,Stim.repeatNum);
Cloud.ScoreMap=zeros(Nx,Ny,Nz);
handles.Setup.DMD.SequenceControl.RepeatModeValue=1;
% ToSave.IndSeq=zeros(Stim.repeatNum,NumofSeq);
Stim.IndSeq=Cloud.IndSeq;
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;

for v=1:numel(Stim.Voltageramp)
    if status == 0
        disp('Procedure interrupted'); break; 
    elseif status ==2
        f = figure;set(f,'Position',[800, 800, 400,200]);
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                                  'Callback','uiresume(gcbf)');
        uiwait(gcf); 
        status = 1;
        close(f);
    end
    if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
       Stim.Output(:,1)=0;
       Stim.Output(:,2) = Stim.Voltageramp(v).*Stim.Array';
       Stim.nonzerovalues = find(Stim.Output(:,2));
    else   
       Stim.Output(:,1) = Stim.Voltageramp(v).*Stim.Array';
       Stim.Output(:,2)=0;
       Stim.nonzerovalues = find(Stim.Output(:,1));
    end
for j = 1:Stim.repeatNum %do as many repetitions as needed 
    if status == 0
        disp('Procedure interrupted'); break; 
    elseif status ==2
        f = figure;set(f,'Position',[800, 800, 400,200]);
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                      'Callback','uiresume(gcbf)');
        uiwait(gcf); 
        status = 1;
        close(f);
    end
    handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
    % generate random combination of masks
%     if v==1
%         Stim.IndSeq(:,j)=randperm(Nx*Ny*Nz);
%         for k=1:NumofSeq
%                 DMDFrames_final = false(handles.Setup.DMD.LX,handles.Setup.DMD.LY,seqL);
%                 for i=1:seqL/handles.Setup.PointCloud.divider
%                     tempInd=Stim.IndSeq((1+(k-1)*seqL*Stim.NumofSpot/handles.Setup.PointCloud.divider+(i-1)*Stim.NumofSpot):((k-1)*seqL*Stim.NumofSpot/handles.Setup.PointCloud.divider+i*Stim.NumofSpot),j);
%                     for t=1:handles.Setup.PointCloud.divider
%                         DMDFrames_final(:,:,t+(i-1)*handles.Setup.PointCloud.divider)=full(reshape(sum(DMDFrames.(sprintf('s%d',t))(:,tempInd),2),handles.Setup.DMD.LX,handles.Setup.DMD.LY));
%                     end
%                 end
%                 [handles.Setup,handles.sequenceid3DXYZ(k+(j-1)*NumofSeq)] = function_StoreImages_DMD(handles.Setup, uint8(DMDFrames_final)*255);
%                 disp(['--loading mask: ' num2str(k) '/' num2str(NumofSeq) ' finished']);
%         end
%     end

    for i=1:NumofSeq
        if status == 0
            disp('Procedure interrupted'); break; 
        elseif status ==2
            f = figure;set(f,'Position',[800, 800, 400,200]);
            h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                          'Callback','uiresume(gcbf)');
            uiwait(gcf); 
            status = 1;
            close(f);
        end
        handles.sequenceid3DXYZ=Cloud.sequenceid3DXYZ;
        Currentsequenceid=handles.sequenceid3DXYZ(i+(j-1)*NumofSeq);
        handles.Setup = function_StartProj_DMD(handles.Setup, Currentsequenceid);
%         handles.Setup = function_feed_DMD(handles.Setup, uint8(DMDFrames_final)*255);
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
        DataSave.XYZ{i,j,v}=Data;
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
        axes(handles.axes5);plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage');title([num2str(seqL) ' spot measurements']);
        disp(['Finish sequence #' num2str(i) '/' num2str(NumofSeq) ' of repeat #' num2str(j) ' at voltage #' num2str(v)]);
    end
%     for i=1:NumofSeq
%         if status == 0
%             disp('Procedure interrupted'); break; 
%         elseif status ==2
%             f = figure;set(f,'Position',[800, 800, 400,200]);
%             h = uicontrol('Position',[20 20 200 40],'String','Continue',...
%                           'Callback','uiresume(gcbf)');
%             uiwait(gcf); 
%             status = 1;
%             close(f);
%         end
%         Data_temp=DataSave.XYZ{i,j}.*Stim.CropMask';
%         Stim.baseline=mean(DataSave.XYZ{i,j}(1:50));
%         Data_temp(Data_temp==0)=[];
%         if rem(numel(Data_temp),seqL)~=0
%            Data_temp=cat(1,Data_temp,0);
%         end
%         DataMatrix=reshape(Data_temp,numel(Data_temp)/seqL,seqL);
%         DataMatrix=DataMatrix';
%         if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
%             for k=1:seqL/handles.Setup.PointCloud.divider
%                 tempInd=Stim.IndSeq((1+(i-1)*seqL*Stim.NumofSpot/handles.Setup.PointCloud.divider+(k-1)*Stim.NumofSpot):((i-1)*seqL*Stim.NumofSpot/handles.Setup.PointCloud.divider+k*Stim.NumofSpot),j);
%                 [indx,indy,indz]=ind2sub([Nx,Ny,Nz],tempInd);
%                 for n=1:Stim.NumofSpot
%                 [ Score4D(indx(n),indy(n),indz(n),j), ~] = function_scorespikes( handles.Setup,Stim,DataMatrix(k,:));
%                 end
%             end
%     end

    disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' at voltage#' num2str(v) ' finished!']);
end
    if status == 0
        disp('Procedure interrupted'); break; 
    elseif status ==2
        f = figure;set(f,'Position',[800, 800, 400,200]);
        h = uicontrol('Position',[20 20 200 40],'String','Continue',...
                      'Callback','uiresume(gcbf)');
        uiwait(gcf); 
        status = 1;
        close(f);
    else
    % plot the result
    yc=zeros(Stim.Npeaks*NumofSeq,Stim.repeatNum);
    for j=1:Stim.repeatNum
        for i=1:NumofSeq
            Data_temp=DataSave.XYZ{i,j,v}.*Stim.CropMask';
            Stim.baseline=mean(DataSave.XYZ{i,j,v}(1:50));
            Data_temp(Data_temp==0)=[];
            if rem(numel(Data_temp),Stim.Npeaks)~=0
               Data_temp=cat(1,Data_temp,0);
            end
            DataMatrix=reshape(Data_temp,numel(Data_temp)/Stim.Npeaks,Stim.Npeaks);
            DataMatrix=DataMatrix';
            for k=1:Stim.Npeaks
                [yc(k+(i-1)*Stim.Npeaks,j),~]=function_scorespikes( handles.Setup,Stim,DataMatrix(k,:));
            end
        end
    end
    h=zeros(Nx*Ny*Nz/Stim.NumofSpot*Stim.repeatNum,Nx*Ny*Nz);
    Ind=Stim.IndSeq(:);
    for i=1:size(h,1)   
        tempInd=Ind((1+(i-1)*Stim.NumofSpot):i*Stim.NumofSpot);
        h(i,tempInd)=1;
    end
    
    %FISTA+TV
    rng('shuffle');
    x0=rand(Nx*Ny*Nz,1);
    h_T=conj(h');
    L=2000;
    alpha=0.005;
    lam=0.5;
    max_iter = 800;
    xk=x0;
    tk=1;%initial tk
    u=xk;
    Yc=yc(:);
        for i=1:max_iter
            Uk=u;
            residue=h*Uk-Yc;
            gradf = h_T*residue;
            Uk1 = Uk - 1/L*gradf;
            Uk13D=reshape(Uk1,[Nx,Ny,Nz]);
            K=2*Uk13D-circshift(Uk13D,1,1)-circshift(Uk13D,1,2);
            xk1=wthresh(Uk1-alpha*K(:),'s',lam/L);
            tk1=1+sqrt(1+4*tk^2)/2;
            u=xk1+(tk-1)/tk1*(xk1-xk);
            xk=xk1;
            if i==1
                loss=sqrt(sum((h*xk-Yc).^2))+lam*sum(abs(xk(:)));
            else
                loss=cat(1,loss,sqrt(sum((h*xk-Yc).^2))+lam*sum(abs(xk(:))));
            end
            if i>2 && loss(end)>loss(end-1)
                break;
            end
        end
    xfinal=reshape(xk,[Nx,Ny,Nz]);  
    
    axes(handles.axes2); imagesc(UY, UX, max(xfinal,[],3)'); xlabel('X pixels'); ylabel('Y pixels'); title(['XYZ map: MZP' num2str(j)]);colorbar;

    figure(800);
    set(gcf, 'Position',  [100, 100, Nz*150, numel(Stim.Voltageramp)*150]);
    for nz=1:Nz
        subplot(numel(Stim.Voltageramp),Nz,(nz+Nz*(v-1)));  
        imagesc(UX,UY,xfinal(:,:,nz)');
        title([' Z=' num2str(UZ(nz)) 'V=' num2str(Stim.Voltageramp(v))]);
        caxis([min(xfinal(:)) max(xfinal(:))]);
        axis image;
        xlabel('x(\mum)');ylabel('y(\mum)');  
    end
    axes(handles.axes3);plot(loss);xlabel('iteration #');ylabel('loss');
    disp(['Voltage #' num2str(v) '/' num2str(numel(Stim.Voltageramp)) ' finished!']);
    end
end
Cloud.ScoreMap=xfinal;
% if isempty(find(ScoreMap, 1))
%    disp('No spike is found in the mapping result!');
% else
%     [counts,edges] = histcounts(ScoreMap);
%     thres=numel(counts)-1;
%     while sum(counts(thres:end))<(0.1*Nx*Ny*Nz)
%         thres=thres-1;
%         if thres<1
%             break;
%         end
%     end
% %     temp=sub2ind(size(ScoreMap),xind,yind, zind);
% %     [Indnonzero,timeSeq,~]=intersect(temp,find(ScoreMap>edges(thres+1)),'stable');
%     [Indnonzero,timeSeq,~]=intersect(1:Nx*Ny*Nz,find(ScoreMap>edges(thres+1)),'stable');
%     timeMap=zeros(size(ScoreMap));
%     timeMap(Indnonzero)=timeSeq;
%     figure(300);set(gcf, 'Position',  [100, 600, 500, 400]);imshow3D(timeMap);title('temporal sequence');colorbar;
% end
% figure(100); set(gcf,'position',[50,50,Nz*150,150]);
% for zz=1:Nz
%    subplot(1,Nz,zz);imagesc(log(ScoreMap(:,:,zz)));axis image;   
% end
% ScoreMap(ScoreMap<=edges(thres+1))=0;%remove background
% Cloud.ScoreMap=ScoreMap;
% figure(101); set(gcf,'position',[50,200,Nz*150,150]);
% for zz=1:Nz
%    subplot(1,Nz,zz);imagesc(Cloud.ScoreMap(:,:,zz));axis image;   
% end
try
ToSave.type = 'Digital random mapping 3D';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
ToSave.Stim.NumofSpot=NumofSpot;
ToSave.Stim.NumofSeq=NumofSeq;
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
catch; 
end;

% ToSave.Score = Score4D;
% figure(200);imshow3D(mean(Score4D,4));
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end

end

function edit56_Callback(hObject, eventdata, handles)
end

% --- # of foci to stimulate simultaneously
function edit56_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Multi-targets illumination simultaneously from previous random
% mapping
function pushbutton52_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'Digital random mapping';
outputSingleScan(handles.Setup.Daq,DAQstate);
% set(handles.slider3,'Value',0);
% set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);

ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

if isempty(find(Cloud.ScoreMap, 1))
    disp('No spike is found in the mapping result!');
else
    ScoreMap=Cloud.ScoreMap;
    [Indx,Indy,Indz]=ind2sub(size(ScoreMap),find(ScoreMap));
    for zz=1:size(ScoreMap,3)
        figure(600);subplot(1,size(ScoreMap,3),zz);
        imagesc(ScoreMap(:,:,zz)');
        axis image;
    end
    set(gcf,'position', [50,550,size(ScoreMap,3)*150,150]);
end

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;

X=UX(Indx);
Y=UY(Indy);
Z=UZ(Indz);
xlswrite('autotargets.xls', [X',Y',Z']);
disp('ROI locations saved in autotargets.xls!');

DataSave.X=X;
DataSave.Y=Y;
DataSave.Z=Z;

Stim.Npeaks = 1;% instead of repeat one holograph, change holograph equals to the sequence length
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array1=Stim.Array;
Stim.CropMask=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array1 = Stim.Array1+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.CropMask = Stim.CropMask+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+i/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox8,'Value') == 1 || get(handles.checkbox11,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.nonzerovalues = find(Stim.Output(:,2));
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
    Stim.nonzerovalues = find(Stim.Output(:,1));
end
Stim.Output(:,5)=0;
select = Stim.UUT<handles.Setup.Scorepikes.sealtestduration;
Stim.Output(select,5) = handles.Setup.Scorepikes.sealtestvalue;
Stim.Output(:,6)=Stim.Baseline.*Stim.Array1';
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array1';
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
[PXX,PYY,PZZ] = ndgrid(UX,UY,UZ);
DataSave.PXX = PXX(:);
DataSave.PYY = PYY(:);
DataSave.PZZ = PZZ(:);

DataSave.XYZ = {}; 
handles.Setup.DMD.SequenceControl.RepeatModeValue=1;

axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
 for j = 1:floor(str2num(get(handles.edit9,'string'))) %do as many repetitions as needed
     if status == 0; disp('Procedure interrupted'); break; end
     if numel(X)>500 disp('Too many ROIs!'); break; end
        [DMDFrames] =  function_makespots_ori(handles.Setup,X,Y,Z,Stim.TargetRadius*ones(size(X)));
        handles.Setup = function_feed_DMD(handles.Setup, uint8(gather(DMDFrames))*255);
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
        Stim.baseline=mean(Data(1:100));
%         if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else handles.Setup.Scorepikes.Method=0; end;
%         [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
%         axes(handles.axes2);  scatter(j, score,'red','filled'); hold on; xlabel('repeat#'); ylabel('Score'); title(['Stimulate ' num2str(indcomb) ' spots']);
        axes(handles.axes5);  plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage');title([num2str(numel(X)) ' spots stim']);
        DataSave.XYZ{j} = Data;
        handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
 end
try
ToSave.type = 'Multi-target illu based on digital mapping';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);
catch 
end
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);

if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end
end


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
end



function edit57_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function edit57_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in 'reshuffle each repetition'.
function checkbox14_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in 'randomize power'.
function checkbox15_Callback(hObject, eventdata, handles)
end


% --- Executes on image panel: preview.
function pushbutton54_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
handles.Setup.src.TriggerMode = 'Internal Trigger';
triggerconfig(handles.Setup.camera, 'immediate');
% handles.Setup.src.AutoContrast = 'ON-MANUAL';
handles.Setup.src.PortSpeedGain = 'Port0-Speed1-100MHz-16bit-Gain1-HDR';
% handles.Setup.src.AutoContrastMax = 3000;
handles.Setup.src.AutoContrast = 'ON-AUTO';
handles.Setup.src.Exposure = str2num(get(handles.edit58,'string'));
ROIsize=round(1200/str2num(get(handles.edit65,'string')));
handles.Setup.camera.ROIPosition = [round((1200-ROIsize)/2) round((1200-ROIsize)/2) ROIsize ROIsize];
preview(handles.Setup.camera);
end


% camera exposure time
function edit58_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function edit58_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- snap images with internal trigger.
function pushbutton55_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global SaveID; SaveID=SaveID+1;
global Cloud;

ToSave.DMDSTATE = 'Snap images';
handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end

Stim.Voltage=str2num(get(handles.edit3,'string'));
Stim.VoltageCa=str2num(get(handles.edit4,'string'));
DAQstate = [Stim.Voltage Stim.VoltageCa handles.Setup.PointCloud.GalvoOffsetVoltage(1) handles.Setup.PointCloud.GalvoOffsetVoltage(2) 0 0 1 0];
outputSingleScan(handles.Setup.Daq,DAQstate);

handles.Setup.src.TriggerMode = 'Internal Trigger';
handles.Setup.src.AutoContrast = 'OFF';
handles.Setup.camera.ROIPosition = [0 0 1200 1200];
triggerconfig(handles.Setup.camera, 'immediate');
handles.Setup.src.Exposure = str2num(get(handles.edit58,'string'));
handles.Setup.camera.FramesPerTrigger = str2num(get(handles.edit60,'string'));

start(handles.Setup.camera);
Images = squeeze(getdata(handles.Setup.camera));
stop(handles.Setup.camera);

Stim.Camera=handles.Setup.camera;
Stim.CameraSrc=handles.Setup.src;
DAQstate = handles.Setup.DAQstateZero;
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
outputSingleScan(handles.Setup.Daq,DAQstate);
figure(255);imshow3D(Images);
ToSave.WF=Images;
ToSave.Stim=Stim;
if get(handles.checkbox1,'Value') == 1
    try
    ToSave.DataSave = {};
    catch
    end
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave','-v7.3');
    disp('Data Saved')
end

end

% --- Executes on button press in pushbutton60: snap a 3D image stack with
% machanical scanning and prime 95B camera
function pushbutton60_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
% DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global SaveID; SaveID=SaveID+1;
global Cloud;

ToSave.DMDSTATE = 'Snap a 3D image stack';
handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end


% stage    
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))))/16;
Stim.CurrentXYZ = getPosition(handles.Setup.SutterStage);
disp(['Current stage Z position:' num2str(Stim.CurrentXYZ(3))]);

%camera
handles.Setup.src.TriggerMode = 'Internal Trigger';
handles.Setup.camera.ROIPosition = [0 0 1200 1200];
triggerconfig(handles.Setup.camera, 'immediate');
handles.Setup.src.Exposure = str2num(get(handles.edit58,'string'));
handles.Setup.src.AutoContrast = 'OFF';
handles.Setup.camera.FramesPerTrigger = 1;
axes(handles.axes3); hold off; cla;

% DMD
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
end

if get(handles.checkbox18,'Value') == 1 
    % generate stripe patterns
    Stim.Grid_period = 20;
    a=rem(1:1200,Stim.Grid_period);
    a=a>Stim.Grid_period/2;
    SIimage=repmat(a,[1200,1]);
    load('AffineMatrix.mat'); % mytform, generated by register button
    DMDsize = imref2d([handles.Setup.DMD.LX,handles.Setup.DMD.LY]);
    SI = uint8(imwarp(SIimage,mytform,'OutputView',DMDsize))*255;
    [fw,fhp,flp]=function_CreatHiLoFilter(SIimage,Stim.Grid_period);
end
HOMO = ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY,'uint8')*255;
ToSave.DataSave={};

Stim.Voltage=str2num(get(handles.edit3,'string'));
Stim.VoltageCa=str2num(get(handles.edit4,'string'));
DAQstate = [Stim.Voltage Stim.VoltageCa handles.Setup.PointCloud.GalvoOffsetVoltage(1) handles.Setup.PointCloud.GalvoOffsetVoltage(2) 0 0 1 0];

function_directfeed_DMD( handles.Setup,HOMO);
outputSingleScan(handles.Setup.Daq,DAQstate);
pause(0.1);
for nz = 1:numel(UZ)
    if status == 0; disp('Procedure interrupted'); break; end;
    moveTime = moveTo(handles.Setup.SutterStage,[Stim.CurrentXYZ(1);Stim.CurrentXYZ(2);Stim.CurrentXYZ(3)+UZ(nz)]);
    pause(0.1);
    start(handles.Setup.camera);
    if nz==1
        HOMOImages = squeeze(getdata(handles.Setup.camera));
    else
        HOMOImages = cat(3,HOMOImages,squeeze(getdata(handles.Setup.camera)));
    end
    stop(handles.Setup.camera);

%     if get(handles.checkbox18,'Value') == 1 
%         function_directfeed_DMD( handles.Setup,SI);
%         outputSingleScan(handles.Setup.Daq,DAQstate);
%         start(handles.Setup.camera);
%         if nz==1
%             SIImages = squeeze(getdata(handles.Setup.camera));
%         else
%             SIImages = cat(3,SIImages,squeeze(getdata(handles.Setup.camera)));
%         end
%         stop(handles.Setup.camera);
% 
%         if nz==1
%             Images = function_HiLo2D(HOMOImages(:,:,end), SIImages(:,:,end), Stim.Grid_period,fw,fhp,flp);
%         else
%             Images = cat(3, Images, function_HiLo2D(HOMOImages(:,:,end), SIImages(:,:,end), Stim.Grid_period,fw,fhp,flp));
%         end  
%     end
end
if get(handles.checkbox18,'Value') == 0
    Images = HOMOImages;
    ToSave.WF = max(Images,[],3);
else
    Images =abs(Images);
    ToSave.WF = max(HOMOImages,[],3);
end
DAQstate = handles.Setup.DAQstateZero;
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
outputSingleScan(handles.Setup.Daq,DAQstate);

figure(255);imshow3D(Images);

axes(handles.axes3); imagesc(ToSave.WF);colorbar;title('MZP');
Stim.UZ=UZ;
Stim.Camera=handles.Setup.camera;
Stim.CameraSrc=handles.Setup.src;
moveTime = moveTo(handles.Setup.SutterStage,[Stim.CurrentXYZ(1);Stim.CurrentXYZ(2);Stim.CurrentXYZ(3)]);
if get(handles.checkbox18,'Value') == 1 
    ToSave.DataSave{1} = Images;
    ToSave.DataSave{2} = HOMOImages;
    ToSave.DataSave{3} = SIImages;
else
    ToSave.DataSave{1} = Images;
end
ToSave.Stim = Stim;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave','-v7.3');
    disp('Data Saved');
end
end

% # of images
function edit60_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function edit60_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- selected ROI simultaneous illumination 2D -----------
function pushbutton56_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = '2D all-optical interrogation';
ToSave.Score=[];
ToSave.Data=[];

% Generate illumination trigger----------------------------------------
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);

Stim.Npeaks = str2num(get(handles.edit6,'string')); % Number of blue light pulses in test
Stim.repeatNum=floor(str2num(get(handles.edit9,'string')));% number of stimulation pulses     
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.VoltageCa = str2num(get(handles.edit63,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),...
    floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
if numel(Stim.Voltageramp)>8 % means forget to change voltage setting back from optimal voltage test
    Stim.Voltageramp=Stim.Voltage;
end
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.PowerMask=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Array=circshift(Stim.Array',-round(Stim.DelayBorder*0.99*handles.Setup.Daq.Rate));
Stim.Baseline = Stim.Output(:,1);
% select = Stim.UUT<handles.Setup.Scorepikes.sealtestduration;
% Stim.Output(select,8) = handles.Setup.Scorepikes.sealtestvalue;
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array;
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));

Stim.Output(:,1)=Stim.Array.*Stim.Voltage;
% Stim.Output(:,1)=0;
Stim.Output(:,8)=Stim.Array.*str2num(get(handles.edit66,'string'));

% Camera settings--------------------------------------------------------
preview='off';
Stim.Exposure = str2num(get(handles.edit58,'string')); %double
if strcmp(preview,'off')
    handles.Setup.src.TriggerMode = 'Trigger first'; % continue capture, Edge Trigger, Trigger first
    % handles.Setup.camera.FramesPerTrigger = Inf;
    triggerconfig(handles.Setup.camera, 'hardware', 'Falling edge', 'Extern');
    handles.Setup.src.Exposure = Stim.Exposure; %int32 type
    handles.Setup.src.AutoContrast = 'OFF';
    handles.Setup.src.FanSpeed = 'Low';
    handles.Setup.src.ExposeOutMode = 'All Rows'; %First Row
    handles.Setup.src.PortSpeedGain = 'Port0-Speed1-100MHz-16bit-Gain1-HDR';
    handles.Setup.src.ClearMode = 'Post-Sequence';
end

if Stim.repeatNum==1 && Stim.FreqHZ<10 % long stim
    handles.Setup.camera.FramesPerTrigger = floor((1000/Stim.FreqHZ-Stim.DurationMS)*0.99/Stim.Exposure);
elseif Stim.repeatNum>1 && Stim.FreqHZ>=10 % pulse train
    handles.Setup.camera.FramesPerTrigger = floor(Stim.DelayBorder*1.8*1000/Stim.Exposure);%600
else
    handles.Setup.camera.FramesPerTrigger = 1;
end
% %edge trigger-------------------------------------------
% % Stim.CropMaskCam=Stim.UUT-Stim.UUT;
% % for i = 1:handles.Setup.camera.FramesPerTrigger
% %     Stim.CropMaskCam = Stim.CropMaskCam+double(Stim.UUT>Stim.DelayBorder+(i-1)*Stim.Exposure/500)...
% %         .*double(Stim.UUT<Stim.DelayBorder+Stim.Exposure/1000+(i-1)*Stim.Exposure/500);
% % end
% % Stim.Output(:,5)=Stim.CropMaskCam;
% trigger first
Stim.CropMask2=zeros(size(Stim.UUT))+double(Stim.UUT>Stim.DelayBorder).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(Stim.Npeaks-1)/Stim.FreqHZ);
Stim.Output(:,5)=circshift(ones(size(Stim.CropMask2'))-Stim.CropMask2',-round(Stim.DelayBorder*0.95*handles.Setup.Daq.Rate)); % trigger camera, 0.4
Stim.Output(1:round((Stim.DurationMS/1000+(Stim.Npeaks-1)/Stim.FreqHZ)*handles.Setup.Daq.Rate),5)=0;

Stim.Output(:,2) = Stim.VoltageCa;
% Stim.Output(:,2)=Stim.Array.*Stim.VoltageCa+3;
% Stim.Output(:,2)=0;

if get(handles.checkbox4,'Value') == 1
    Stim.Output(:,3)=handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1);
    Stim.Output(:,4)=handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2);
%     Stim.Output(:,3)=handles.Setup.PointCloud.GalvoOffsetVoltageBlue(1);
%     Stim.Output(:,4)=handles.Setup.PointCloud.GalvoOffsetVoltageBlue(2);
else
    Stim.galvomask=circshift(ones(size(Stim.CropMask2'))-Stim.CropMask2',-round(Stim.DelayBorder*0.99*handles.Setup.Daq.Rate));
    Stim.Output(:,3)=Stim.Output(:,3).*Stim.Array+...
        Stim.galvomask.*handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1);
    Stim.Output(:,4)=Stim.Output(:,4).*Stim.Array+...
        Stim.galvomask.*handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2);
end

[idy,idx] = find(ToSave.ImageMask);
Stim.Padarray=0;
Stim.ROIoffset=[min(idx)-Stim.Padarray min(idy)-Stim.Padarray max(idx)-min(idx)+Stim.Padarray max(idy)-min(idy)+Stim.Padarray];
handles.Setup.camera.ROIPosition = Stim.ROIoffset;

% stage    
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))))/16;
Stim.CurrentXYZ = getPosition(handles.Setup.SutterStage);
disp(['Current stage Z position:' num2str(Stim.CurrentXYZ(3))]);
Nz=numel(UZ);
if Nz>10
    UZ=0;
    Nz=1;
end

% calculate ROI pattern on DMD coordinate------------------------------------------
load('AffineMatrix.mat'); % mytform, generated by register button
DMDsize = imref2d([handles.Setup.DMD.LX,handles.Setup.DMD.LY]);
Stim.ROI = imwarp(ToSave.ImageMask,mytform,'OutputView',DMDsize);
Stim.StimROI = imwarp(ToSave.StimMask,mytform,'OutputView',DMDsize);
axes(handles.axes2); hold off; cla;axes(handles.axes5); hold off; cla;
DataSave.Images={};DataSave.Iref={};DataSave.UZ=UZ*16;
Mask=ToSave.ImageMaskLabel(Stim.ROIoffset(2):Stim.ROIoffset(2)+Stim.ROIoffset(4)-1, ...
     Stim.ROIoffset(1):Stim.ROIoffset(1)+Stim.ROIoffset(3)-1);
Stim.ROILabel=bwlabel(Stim.StimROI, 8);
% N=max(Stim.ROILabel(:));
N=max(Mask(:));
for ii=1:N
    ROIarea(ii)=numel(nonzeros(Mask==ii));
end
Cloud.permutepattern=[];

if Stim.repeatNum==1 && Stim.FreqHZ<10 % long stim
    Stim.cutoffFreq=Stim.FreqHZ;
elseif Stim.repeatNum>1 && Stim.FreqHZ>=10 % pulse train
    Stim.cutoffFreq=1/(Stim.DelayBorder*2);
else
    Stim.cutoffFreq=0.5; % other condition
end

% Project patterns
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
end
function_directfeed_DMD( handles.Setup,zeros(handles.Setup.DMD.LX,handles.Setup.DMD.LY,'uint8'));
outputSingleScan(handles.Setup.Daq,Stim.Output(1,:));
pause(30); % wait the yellow laser to be stable;
if get(handles.checkbox16,'Value') == 0
    Stim.permute='false';
    % illuminate all ROIs together
    Mask3D=repmat(logical(ToSave.ImageMask(Stim.ROIoffset(2):Stim.ROIoffset(2)+Stim.ROIoffset(4)-1, ...
     Stim.ROIoffset(1):Stim.ROIoffset(1)+Stim.ROIoffset(3)-1)),[1,1,handles.Setup.camera.FramesPerTrigger*Stim.repeatNum,Nz]);
    Images_all=zeros(size(Mask,1),size(Mask,2),handles.Setup.camera.FramesPerTrigger*Stim.repeatNum,Nz,'uint16');
    function_directfeed_DMD( handles.Setup,uint8(gather(Stim.ROI(:,:,1)))*255);
    for nz = 1:Nz
        if status == 0; disp('Procedure interrupted'); break; end;
        moveTime = moveTo(handles.Setup.SutterStage,[Stim.CurrentXYZ(1);Stim.CurrentXYZ(2);Stim.CurrentXYZ(3)+UZ(nz)]);
        pause(0.1);
    % Ready to collect calcium data --------------------------------------
%         start(handles.Setup.camera);
%         for j=1:Stim.repeatNum
%             if status == 0; disp('Procedure interrupted'); break; end;
%             queueOutputData(handles.Setup.Daq,Stim.Output);
%             DataSave.Iref{j}=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
%             axes(handles.axes5);plot(Stim.UUT, medfilt1(DataSave.Iref{j}, 500));title('Ref: laser intensity');xlabel('time (s)');
%         end
%         Images=squeeze(getdata(handles.Setup.camera));
%         stop(handles.Setup.camera);
   
    %-trigger readout---------------------------------------------------
        for j=1:Stim.repeatNum
            start(handles.Setup.camera);
            if status == 0; disp('Procedure interrupted'); break; end;
            queueOutputData(handles.Setup.Daq,Stim.Output);
            Iref=startForeground(handles.Setup.Daq);
            if j==1
                %RawData is the analog input signal from AI16 pin
%                 Iref1=nonzeros(Iref.*Stim.Output(:,5));
                Iref1=Iref;
            else
%                 Iref1=cat(1,Iref1,nonzeros(Iref.*Stim.Output(:,5)));
                Iref1=cat(1,Iref1,Iref);
            end
            T_ref=0:1/handles.Setup.Daq.Rate*1000:(numel(Iref1)-1)/handles.Setup.Daq.Rate*1000;
            axes(handles.axes5);plot(T_ref,Iref1);xlabel('time [ms]');title('Ref: yellow laser intensity');
            stop(handles.Setup.camera);
            Images_all(:,:,(j-1)*handles.Setup.camera.FramesPerTrigger+1:j*handles.Setup.camera.FramesPerTrigger,nz)=squeeze(getdata(handles.Setup.camera));
            DataSave.Iref{j,nz}=Iref;
            disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
        end
        disp(['nz=' num2str(nz) '/' num2str(Nz) ' finished!']);
    end
    moveTime = moveTo(handles.Setup.SutterStage,[Stim.CurrentXYZ(1);Stim.CurrentXYZ(2);Stim.CurrentXYZ(3)]);
    DAQstate = handles.Setup.DAQstateZero;
    outputSingleScan(handles.Setup.Daq,DAQstate);
    %processing data------------------------------------------------------
    disp('start processing data...');
    %calcium trace continue recording--------------------------------
%     Ica=zeros(handles.Setup.camera.FramesPerTrigger*Stim.repeatNum, Nz, N);
%     figure(500);set(gcf,'position',[100,100,ceil(N/10)*200,100*10]);
%     for ii=1:N
%         masktemp=Mask3D==ii;
%         c=nonzeros(double(Images_all).*double(masktemp));
%         DataSave.Images{ii}=reshape(c,[numel(c)/(handles.Setup.camera.FramesPerTrigger*Stim.repeatNum*Nz),...
%         handles.Setup.camera.FramesPerTrigger*Stim.repeatNum,Nz]);
%         b=squeeze(mean(DataSave.Images{ii},1));
%         Ica(:,:,ii)=b;
%         figure(500);subplot(10,ceil(N/10),ii);plot(b);axis tight;
%         drawnow;
%     end
%--------------------------------------------------------------
    c=nonzeros(double(Images_all).*double(Mask3D));
    Non0=reshape(c,[numel(c)/(handles.Setup.camera.FramesPerTrigger*Stim.repeatNum*Nz),...
        handles.Setup.camera.FramesPerTrigger*Stim.repeatNum,Nz]);
    Ica=zeros(handles.Setup.camera.FramesPerTrigger*Stim.repeatNum, Nz, N);
    figure(500);set(gcf,'position',[100,100,ceil(N/10)*200,100*10]);
    for ii=1:N
        if ii==1
            DataSave.Images{ii}=Non0(1:ROIarea(1),:,:);
        else
            DataSave.Images{ii}=Non0((1+sum(ROIarea(1:ii-1))):sum(ROIarea(1:ii)),:,:);
        end
        b=squeeze(mean(DataSave.Images{ii},1));
        Ica(:,:,ii)=b;
        figure(500);subplot(10,ceil(N/10),ii);plot(b);axis tight;
        drawnow;
    end
    %compute dF/F of the brightest ROI-----------------------------
    [~,indmaxlinear]=max(Ica(:));
    [~,idy,indmax]=ind2sub(size(Ica),indmaxlinear);
    CaTrace=Ica(:,idy,indmax)-100;
    dF=zeros(size(CaTrace));
    w=round((1000/Stim.Exposure)/Stim.cutoffFreq*0.5);
    for ii=1:numel(CaTrace)
        if ii<=w
            F0=min(CaTrace(1:ii+w));
        elseif ii>numel(CaTrace)-w
            F0=min(CaTrace(ii-w:end));
        else
            F0=min(CaTrace(ii-w:ii+w));
        end
        dF(ii)=(CaTrace(ii)-F0)/F0;
    end
    T_ds=0:Stim.Exposure:(numel(CaTrace)-1)*Stim.Exposure;
    axes(handles.axes2);plot(T_ds,CaTrace);title(['Ca2+ Intensity of Imax ROI#' num2str(indmax) ', Z=' num2str(DataSave.UZ(idy))]);xlabel('time [ms]');ylabel('Intensity');
    axes(handles.axes5);plot(T_ds,dF);xlim([min(T_ds) max(T_ds)]);xlabel('time [ms]');ylabel('\DeltaF/F');
    title(['\DeltaF/F of of Imax ROI#' num2str(indmax) ', Z=' num2str(DataSave.UZ(idy)) '\mum']);
    figure(600);imshow3D(Images_all(:,:,:,idy));
else
    % permute
    Stim.permute='true';
    temp = double(Stim.UUT>Stim.DelayBorder).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(Stim.Npeaks-1)/Stim.FreqHZ);
    Stim.Output(:,6)=circshift(temp',-round(Stim.DelayBorder*0.99*handles.Setup.Daq.Rate))+Stim.Output(:,5);
    ProjMode = function_CheckProjMode_DMD(handles.Setup);
    if ProjMode==2301
        [handles.Setup]=function_StopProj_DMD(handles.Setup);
        [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
    end
    Stim.Nregions=str2num(get(handles.edit64,'string'));% number of big regions
    DMDFrameROI=zeros(handles.Setup.DMD.LX,handles.Setup.DMD.LY,Stim.Nregions,'uint8');
    %decompose ROI by its label
%     for indcomb=1:Stim.Nregions
%         temp=Stim.ROILabel;
%         temp(temp<=(indcomb-1)* floor(N/Stim.Nregions))=0;
%         if indcomb<Stim.Nregions %last frame includes all ROIs
%             temp(temp>indcomb* floor(N/Stim.Nregions))=0;
%         end
%         temp(temp~=0)=255;
%         DMDFrameROI(:,:,indcomb)=uint8(temp);
%     end
%permute
%     Stim.permROI=randperm(N);
    [~, Stim.permROI] = function_poissonDisc_fixgrid(ToSave.ROIcenter,80);
    Nsubregions=floor(N/Stim.Nregions);
    for indcomb=1:Stim.Nregions
        if indcomb<Stim.Nregions
            indtemp=Stim.permROI(1+(indcomb-1)*Nsubregions:indcomb*Nsubregions);
        else
            indtemp=Stim.permROI(1+(indcomb-1)*Nsubregions:end);
        end
        for kk=1:length(indtemp)
            temp=Stim.ROILabel;
            temp(temp~=indtemp(kk))=0;
            temp(temp~=0)=255;
            DMDFrameROI(:,:,indcomb)=DMDFrameROI(:,:,indcomb)+uint8(temp);
        end
        figure(800);subplot(1,Stim.Nregions,indcomb);imagesc(DMDFrameROI(:,:,indcomb));
        title(num2str(indtemp));axis image;
    end
       
    % compute DMD patterns of all combo
    illuminatetime=Stim.DurationMS*1000+(Stim.Npeaks-1)/Stim.FreqHZ*10^6;%us
%     i=1;
%     for indcomb=1:Stim.Nregions
%         if status == 0; disp('Procedure interrupted'); break; end;
%         Stim.IndComb{indcomb}=nchoosek(1:Stim.Nregions,indcomb);
%         for p=1:size(Stim.IndComb{indcomb},1)
%             ic=Stim.IndComb{indcomb}(p,:);
%             DMDFrames=max(DMDFrameROI(:,:,ic),[],3);
%             [handles.Setup,Cloud.permutepattern(i)] = function_StoreImages_DMD_timecontrol(handles.Setup, DMDFrames, illuminatetime, 1);
%             i=i+1;
%         end
%     end
% [handles.Setup,Cloud.allROIid] = function_StoreImages_DMD_timecontrol(handles.Setup, DMDFrames,...
%     handles.Setup.camera.FramesPerTrigger*Stim.Exposure*1.5*10^3, 1);
    for indcomb=1:Stim.Nregions
        [handles.Setup,Cloud.permutepattern(indcomb)] = function_StoreImages_DMD_timecontrol(handles.Setup, DMDFrameROI(:,:,indcomb), illuminatetime, 1);
    end
    [handles.Setup,Cloud.allROIid] = function_StoreImages_DMD_timecontrol(handles.Setup, uint8(Stim.ROI)*255,...
    handles.Setup.camera.FramesPerTrigger*Stim.Exposure*1.5*10^3, 1);

    Npatterns=numel(Cloud.permutepattern);
    Stim.OutputStim=Stim.Output(1:round(Stim.Npeaks/Stim.FreqHZ*handles.Setup.Daq.Rate),:);
    Stim.OutputImage=Stim.Output(round(Stim.Npeaks/Stim.FreqHZ*handles.Setup.Daq.Rate)+1:end,:);
    for i = 1:Npatterns
        if status == 0; disp('Procedure interrupted'); break; end;
        for j=1:Stim.repeatNum
            if status == 0; disp('Procedure interrupted'); break; end;
            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.permutepattern(i));
            queueOutputData(handles.Setup.Daq,Stim.OutputStim);
            startForeground(handles.Setup.Daq);

            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.allROIid);
            if strcmp(preview,'off')
                start(handles.Setup.camera);
            end
            outputSingleScan(handles.Setup.Daq,[0,Stim.VoltageCa,handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1),...
             handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2), 1, 1, 0, 0]);
            pause(handles.Setup.camera.FramesPerTrigger*Stim.Exposure*1.5/1000);
            if strcmp(preview,'off')
                stop(handles.Setup.camera);
                currentstack=squeeze(getdata(handles.Setup.camera));
                if j==1
                    Images=currentstack(:,:,2:end-2);
                else
                    Images=cat(3,Images,currentstack(:,:,2:end-2));
                end
            end
            disp(['finish ' num2str(j) '/' get(handles.edit9,'string') ' repetition']);
        end
        if strcmp(preview,'off')
            if i==1
                Images_all=Images;
            else
                Images_all=cat(4,Images_all,Images);
            end
        end
        disp(['--------------Pattern#=' num2str(i) '/' num2str(Npatterns) ' finished!']);
    end
    
    handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
    DAQstate = handles.Setup.DAQstateZero;
    outputSingleScan(handles.Setup.Daq,DAQstate);
    %processing data------------------------------------------------------
    disp('start processing data...');
    Nframe=size(Images,3);
    Mask3D=repmat(logical(ToSave.ImageMask(Stim.ROIoffset(2):Stim.ROIoffset(2)+Stim.ROIoffset(4)-1, ...
        Stim.ROIoffset(1):Stim.ROIoffset(1)+Stim.ROIoffset(3)-1)),[1,1,Nframe,Npatterns]);
    c=nonzeros(double(Images_all).*double(Mask3D));
    Non0=reshape(c,[numel(c)/(Nframe*Npatterns),Nframe,Npatterns]);
    Ica=zeros(Nframe, Npatterns, N);
    figure(500);set(gcf,'position',[100,100,ceil(N/10)*200,100*10]);
    for ii=1:N
        if ii==1
            DataSave.Images{ii}=Non0(1:ROIarea(1),:,:);
        else
            DataSave.Images{ii}=Non0((1+sum(ROIarea(1:ii-1))):sum(ROIarea(1:ii)),:,:);
        end
        b=squeeze(mean(DataSave.Images{ii},1));
        Ica(:,:,ii)=b;
        figure(500);subplot(10,ceil(N/10),ii);plot(b(size(b,1)/2:end,:));axis tight;
        drawnow;
    end
    %compute dF/F of the brightest ROI-----------------------------
    [~,indmaxlinear]=max(Ica(:));
    [~,idy,indmax]=ind2sub(size(Ica),indmaxlinear);
    CaTrace=Ica(:,idy,indmax)-100;
    dF=zeros(size(CaTrace));
    w=round((1000/Stim.Exposure)/Stim.cutoffFreq)*0.5;
    for ii=1:numel(CaTrace)
        if ii<=w
            F0=min(CaTrace(1:ii+w));
        elseif ii>numel(CaTrace)-w
            F0=min(CaTrace(ii-w:end));
        else
            F0=min(CaTrace(ii-w:ii+w));
        end
        dF(ii)=(CaTrace(ii)-F0)/F0;
    end
    T_ds=0:Stim.Exposure:(numel(CaTrace)-1)*Stim.Exposure;
    axes(handles.axes2);plot(T_ds,CaTrace);title(['Ca2+ Intensity of Imax ROI#' num2str(indmax) ', pattern#=' num2str(idy)]);xlabel('time [ms]');ylabel('Intensity');
    axes(handles.axes5);plot(T_ds,dF);xlim([min(T_ds) max(T_ds)]);xlabel('time [ms]');ylabel('\DeltaF/F');
    title(['\DeltaF/F of of Imax ROI#' num2str(indmax) ', pattern#=' num2str(idy)]);
    ToSave.DMDFrameROI=DMDFrameROI;
end
disp('Finish data processing!');
Stim.Camera=handles.Setup.camera;
Stim.CameraSrc=handles.Setup.src;
ToSave.DataSave=DataSave;
ToSave.DMDROI=Stim.ROI;
ToSave.ProcessedTraces=Ica;
ToSave.Stim = Stim;
ToSave.AffineMatrix=mytform;
ToSave.status = status;
handles.Setup.src.FanSpeed = 'High';
if get(handles.checkbox1,'Value') == 1
    disp('Start saving...');
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave','-v7.3');
    disp('Data Saved')
end
end

% --- select ROIs.
function pushbutton57_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global SaveID; SaveID=SaveID+1;

ToSave.DMDSTATE = 'Select ROIs';
ToSave.StimMask=[];
ToSave.ImageMask=[];
ToSave.ImageMaskLabel=[];

if length(size(ToSave.WF)) == 3
    I = mean(ToSave.WF,3);
else
    I = ToSave.WF;
end
% if min(I(:))<5 %HiLo
%     Iraw=I;
%     I=zeros(size(Iraw));
%     I(200:1000,200:1000)=Iraw(200:1000,200:1000);
% end
mask1=imbinarize(I,'global');
mask2=imbinarize(I,'adaptive', 'Sensitivity', str2num(get(handles.edit61,'string')));
se = strel('disk',5);
tempMask=imopen(mask1.*mask2,se);
axes(handles.axes4); hold off; cla;
if nnz(tempMask(:))<10*10
    disp('ROI is too small! Please select again!');
else
    tempMaskLabel=bwlabel(tempMask, 8);
    N=max(tempMaskLabel(:));
    %find the center of each ROI
    ToSave.TargetRadius = str2num(get(handles.edit19,'string'));
    [ii,jj] = ndgrid(1:size(I,1),1:size(I,2));
    ToSave.ROIcenter=zeros(N,2);
    StimMask=zeros(size(I));
    for i=1:N
        amask=tempMaskLabel==i;
        tot_mass = sum(amask(:));
        ToSave.ROIcenter(i,1) = round(sum(ii(:).*amask(:))/tot_mass);
        ToSave.ROIcenter(i,2) = round(sum(jj(:).*amask(:))/tot_mass);
        StimMask=StimMask+double(((ii- ToSave.ROIcenter(i,1)).^2+(jj-ToSave.ROIcenter(i,2)).^2)<ToSave.TargetRadius^2);
    end
    StimMask(StimMask>1)=1;%merge overlap regions
    ToSave.StimMask=StimMask;
    ToSave.ImageMask=max(cat(3,StimMask,tempMask),[],3);
    ToSave.ImageMaskLabel=bwlabel(ToSave.ImageMask, 8);
    axes(handles.axes4);
    imagesc(ToSave.ImageMask);title(['Selected ' num2str(N) ' ROIs']);axis image;
    green = cat(3, zeros(size(ToSave.WF)),ones(size(ToSave.WF)),zeros(size(ToSave.WF)));
    hold on
    h=imshow(green);
    hold off
    set(h,'AlphaData',ToSave.StimMask*0.6);
end

if get(handles.checkbox1,'Value') == 1
    try
        ToSave.DataSave = {};
    catch
    end
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end
end


function edit61_Callback(hObject, eventdata, handles)
end

% --- threshold for adaptive image binarization.
function edit61_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- % register DMD with Camera using control points
function pushbutton58_Callback(hObject, eventdata, handles)
global ToSave;
global DAQstate;
global Cloud;

ToSave.DMDSTATE = 'Registration';
% set galvo-mirror to zero and read laser intensity from the control bar
handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
if get(handles.checkbox8,'Value') == 1 %yellow stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageYellow;
elseif get(handles.checkbox11,'Value') == 1 %red stim
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageRed;
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
    Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
end
DAQstate = [DAQstate(1) DAQstate(2) handles.Setup.PointCloud.GalvoOffsetVoltage(1) handles.Setup.PointCloud.GalvoOffsetVoltage(2) 0 0 1 0];

% generate grid pattern
[x,y]=meshgrid(linspace(-300,400,11),linspace(-300,400,11));
DefocusingRadius=zeros(numel(x),1);
TargetRadius=4*abs(x(:).*y(:))/10^5+2;
FinalFrames=function_makespots_ori(handles.Setup,x(:),y(:),DefocusingRadius,TargetRadius);
RegPattern = uint8(gather(FinalFrames(:,:,1))).*255;

ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
end
function_directfeed_DMD( handles.Setup,RegPattern);
outputSingleScan(handles.Setup.Daq,DAQstate);
ToSave.RegPattern=RegPattern;

handles.Setup.src.TriggerMode = 'Internal Trigger';
triggerconfig(handles.Setup.camera, 'immediate');
handles.Setup.camera.ROIPosition = [0 0 1200 1200];
handles.Setup.src.Exposure = 100;
handles.Setup.camera.FramesPerTrigger = 1;
start(handles.Setup.camera);
RegImages = squeeze(getdata(handles.Setup.camera));
ToSave.RegImage=RegImages;
[movingPoints, fixedPoints]=cpselect(RegImages,RegPattern,'wait',true);%remember to change back
movingPointsChanged=cpcorr(movingPoints,fixedPoints,RegImages,RegPattern);%remember to change back
mytform = fitgeotrans(movingPointsChanged, fixedPoints, 'affine');
ToSave.AffineMatrix=mytform;
Stim.Voltage=DAQstate(1);
Stim.VoltageCa=DAQstate(2);
Stim.Camera=handles.Setup.camera;
Stim.CameraSrc=handles.Setup.src;
ToSave.Stim=Stim;
save('AffineMatrix.mat','mytform');
disp('New affine matrix saved!');
% if get(handles.checkbox1,'Value') == 1
%     filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
%     save(filename,'ToSave');
%     disp('Data Saved')
% end
stop(handles.Setup.camera);
end



function edit63_Callback(hObject, eventdata, handles)
end

% --- Ca2+ laser voltage, yellow laser voltage used only for calcium
% imaging
function edit63_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Permute ROIs, =1 all combination of illumination; =0, illuminate all
function checkbox16_Callback(hObject, eventdata, handles)
end


function edit64_Callback(hObject, eventdata, handles)
end

% ---Number of ROIs. Because permute nchoosek can results in a big value,
% this value limits the number of n.
function edit64_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes during object creation, after setting all properties.
function pushbutton55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
end


% --- Executes on button press in checkbox17: readout PPSF by camera.
function checkbox17_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton59: current injection (2.5V) and imaging.
% also can be used for visual stim, just need to change the voltage to 3.3V
function pushbutton59_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'current injection';
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
Stim.Npeaks = floor(str2num(get(handles.edit6,'string'))); % Number of blue light pulses in test
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
Stim.VoltageProbe = str2num(get(handles.edit66,'string')); % voltage for current injection, not for stimulation
Stim.Voltage = str2num(get(handles.edit11,'string')); % blue laser voltage
Stim.VoltageCa = str2num(get(handles.edit63,'string')); % voltage for imaging laser
Stim.TargetRadius = str2num(get(handles.edit19,'string'));
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,7) =double(Stim.Subclock<Stim.DutyCycle);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
Stim.Array2=Stim.Array;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
    Stim.Array2 = Stim.Array2+double(Stim.UUT>Stim.DelayBorder*0.98+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Output(:,1)=Stim.Voltage;
Stim.Output(:,2)=0;
Stim.Output(:,3)=handles.Setup.PointCloud.GalvoOffsetVoltageBlue(1);
Stim.Output(:,4)=handles.Setup.PointCloud.GalvoOffsetVoltageBlue(2);
Stim.Output(:,5)=1;
Stim.Output(:,6)=0;
Stim.Output(:,8) = Stim.VoltageProbe.*Stim.Array';
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;

handles.Setup.src.TriggerMode = 'Trigger first';
triggerconfig(handles.Setup.camera, 'hardware', 'Falling edge', 'Extern');
Stim.Exposure = str2num(get(handles.edit58,'string'));
handles.Setup.src.Exposure = Stim.Exposure;
handles.Setup.src.AutoContrast = 'OFF';
handles.Setup.src.FanSpeed = 'Low';
handles.Setup.src.ExposeOutMode = 'First Row';
handles.Setup.src.PortSpeedGain = 'Port0-Speed1-100MHz-16bit-Gain1-HDR';
handles.Setup.src.ClearMode = 'Post-Sequence';
Stim.ROIoffset=[0 0 1200 1200];% assume the cell is in the center
handles.Setup.camera.ROIPosition = Stim.ROIoffset;
handles.Setup.camera.FramesPerTrigger = Inf;
temp=function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),0,Stim.TargetRadius);
Stim.DMDimage=uint8(gather(temp(:,:,1)))*255;
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
end
function_directfeed_DMD(handles.Setup, Stim.DMDimage);
tempstatus=Stim.Output(1,:);
tempstatus(5)=0;
outputSingleScan(handles.Setup.Daq,tempstatus);
pause(30); % wait the yellow laser to be stable;

start(handles.Setup.camera);
for j = 1:floor(str2num(get(handles.edit9,'string'))); %do as many repetitions as needed
    if status == 0; disp('Procedure interrupted'); break; end;
    queueOutputData(handles.Setup.Daq,Stim.Output);
    if j==1
        Iref=startForeground(handles.Setup.Daq);%RawData is the analog input signal from AI16 pin
    else
        Iref=cat(1,Iref,startForeground(handles.Setup.Daq));
    end
    axes(handles.axes5);plot(medfilt1(Iref, 500));title('ephy read');
end
stop(handles.Setup.camera);
DAQstate = handles.Setup.DAQstateZero;
outputSingleScan(handles.Setup.Daq,DAQstate);

Images=squeeze(getdata(handles.Setup.camera));
mask=zeros(size(Images));
mask(600-Stim.TargetRadius:600+Stim.TargetRadius,600-Stim.TargetRadius:600+Stim.TargetRadius,:)=1;

Non0=nonzeros(double(Images).*mask);
MaskedImages=reshape(Non0,[numel(Non0)/size(Images,3), size(Images,3)]);
Ica=mean(MaskedImages,1);
DataSave.Images=MaskedImages;
DataSave.Iref=Iref;
T_ds=0:Stim.Exposure:(numel(Ica)-1)*Stim.Exposure;
axes(handles.axes2);plot(T_ds,Ica);xlabel('Time [ms]'); ylabel('Intensity');title('Raw Intensty');
CaTrace=Ica-100;%100: background intensity
dF=zeros(size(CaTrace));
Stim.cutoffFreq=1/Stim.DelayBorder;
w=round((1000/Stim.Exposure)/Stim.cutoffFreq);
for ii=1:numel(CaTrace)
    if ii<=w
        F0=min(CaTrace(1:ii+w));
    elseif ii>numel(CaTrace)-w
        F0=min(CaTrace(ii-w:end));
    else
        F0=min(CaTrace(ii-w:ii+w));
    end
    dF(ii)=(CaTrace(ii)-F0)/F0;
end
axes(handles.axes3);plot(T_ds,dF);xlabel('Time [ms]'); ylabel('\DeltaF/F');title('\DeltaF/F');

Stim.Camera=handles.Setup.camera;
Stim.CameraSrc=handles.Setup.src;
figure(258);imshow3D(Images);
ToSave.DataSave=DataSave;
ToSave.ImageMask=mask;
ToSave.ProcessedTraces=Ica;
ToSave.Stim = Stim;
ToSave.status = status;
handles.Setup.src.FanSpeed = 'High';
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave','-v7.3');
    disp('Data Saved')
end  
end



% --- Executes on button press in checkbox18: HiLo de-scattering.
function checkbox18_Callback(hObject, eventdata, handles)
end



function edit65_Callback(hObject, eventdata, handles)
end

% --- Executes zoom-in times of preview
function edit65_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function edit66_Callback(hObject, eventdata, handles)
end

% --- voltage for current injection.
function edit66_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton61: digital PPSF readout with
% camera
function pushbutton61_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'DMD Digital PPSF Sequence read by camera';
% outputSingleScan(handles.Setup.Daq,DAQstate);
% set(handles.slider3,'Value',0);
% set(handles.slider2,'Value',0);
% set(handles.edit3,'string',0);
% set(handles.edit4,'string',0);
Stim.Npeaks = floor(str2num(get(handles.edit6,'string'))); % Number of blue light pulses in test
Stim.repeatNum=floor(str2num(get(handles.edit9,'string')));% number of stimulation pulses     
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/Stim.DurationMS)*Stim.DurationMS);
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/Cloud.CycleLength)*Cloud.CycleLength; % to avoid start from the middle of a trigger
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,7) =0;
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( -str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( -str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( -str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(str2num(get(handles.edit20,'string')),str2num(get(handles.edit47,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(str2num(get(handles.edit21,'string')),str2num(get(handles.edit48,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(str2num(get(handles.edit22,'string')),str2num(get(handles.edit49,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
Nx=floor(str2num(get(handles.edit23,'string')));
Ny=floor(str2num(get(handles.edit24,'string')));
Nz=floor(str2num(get(handles.edit25,'string')));
Cloud.sequenceid1D.x=zeros(Nx,1);
Cloud.sequenceid1D.y=zeros(Ny,1);
Cloud.sequenceid1D.z=zeros(Nz,1);

handles.Setup.src.TriggerMode = 'Trigger first'; % continue capture, Edge Trigger, Trigger first
% handles.Setup.camera.FramesPerTrigger = Inf;
triggerconfig(handles.Setup.camera, 'hardware', 'Falling edge', 'Extern');
Stim.Exposure = str2num(get(handles.edit58,'string')); %double
handles.Setup.src.Exposure = Stim.Exposure; %int32 type
handles.Setup.src.AutoContrast = 'OFF';
handles.Setup.src.FanSpeed = 'Low';
handles.Setup.src.ExposeOutMode = 'All Rows'; %First Row
handles.Setup.src.PortSpeedGain = 'Port0-Speed1-100MHz-16bit-Gain1-HDR';
handles.Setup.src.ClearMode = 'Post-Sequence';
handles.Setup.camera.FramesPerTrigger = floor(Stim.DelayBorder*1.8*1000/Stim.Exposure);%600

Stim.ROIoffset=[500 500 200 200];% assume the cell is in the center
Stim.DMDFrame_forimage=uint8(gather(function_makespots(handles.Setup,0,0,0,100)))*255;
handles.Setup.camera.ROIPosition = Stim.ROIoffset;
Stim.VoltageCa = str2num(get(handles.edit63,'string'));
Stim.Array=circshift(Stim.Array',-round(Stim.DelayBorder*0.99*handles.Setup.Daq.Rate));
Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array;
Stim.VoltageCa = str2num(get(handles.edit63,'string'));
Stim.Output(:,2) = Stim.VoltageCa;
%     Stim.CropMask2=zeros(size(Stim.UUT))+double(Stim.UUT>Stim.DelayBorder+Stim.DurationMS/1000+(Stim.Npeaks-1)/Stim.FreqHZ);
%     Stim.Output(:,5)=Stim.CropMask2';
Stim.CropMask2=zeros(size(Stim.UUT))+double(Stim.UUT>Stim.DelayBorder).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(Stim.Npeaks-1)/Stim.FreqHZ);
Stim.Output(:,5)=circshift(ones(size(Stim.CropMask2'))-Stim.CropMask2',-round(Stim.DelayBorder*0.95*handles.Setup.Daq.Rate)); % trigger camera
Stim.Output(1:round((Stim.Npeaks-1)/Stim.FreqHZ*handles.Setup.Daq.Rate),5)=0;

%     Stim.CropMask2=zeros(size(Stim.UUT))+double(Stim.UUT>Stim.DelayBorder).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+Stim.Npeaks/Stim.FreqHZ);
%     % Stim.Output(:,2)=(ones(size(Stim.Array'))-Stim.Array').*Stim.CropMask2'.*Stim.VoltageCa;
%     Stim.Output(:,2)=Stim.CropMask2'.*Stim.VoltageCa;
Stim.galvomask=circshift(ones(size(Stim.CropMask2'))-Stim.CropMask2',-round(Stim.DelayBorder*0.99*handles.Setup.Daq.Rate));
Stim.Output(:,3)=Stim.Output(:,3).*Stim.Array+...
    Stim.galvomask.*handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1);
Stim.Output(:,4)=Stim.Output(:,4).*Stim.Array+...
    Stim.galvomask.*handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2);
Stim.Output(:,6)=Stim.Baseline.*Stim.Array+Stim.Output(:,5);
Stim.cutoffFreq=1/(Stim.DelayBorder*2);
DataSave.X={};DataSave.SX={};DataSave.Y={};DataSave.SY={};DataSave.Z={};DataSave.SZ={};
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

for i=1:Nx
    if status == 0; disp('Procedure interrupted'); break; end;
    [DMDFrames] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string'))+UX(i),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
    DMDFrames = 255*uint8(gather(DMDFrames));
%             DMDFrames = 255*cat(3,DMDFrames,ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY, 'uint8'));
    illuminatetime=handles.Setup.PointCloud.CycleLength/handles.Setup.PointCloud.divider*10^6*0.6;%us
    [handles.Setup,Cloud.sequenceid1D.x(i)] = function_StoreImages_DMD_timecontrol(handles.Setup, DMDFrames, illuminatetime, Stim.Npeaks);
%             [handles.Setup,Cloud.sequenceid1D.x(i)] = function_StoreImages_DMD(handles.Setup,DMDFrames);
end
for i=1:Ny
    if status == 0; disp('Procedure interrupted'); break; end;
    [DMDFrames] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string'))++UY(i),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
    DMDFrames = 255*uint8(gather(DMDFrames));
%             DMDFrames = 255*cat(3,DMDFrames,ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY, 'uint8'));
    illuminatetime=handles.Setup.PointCloud.CycleLength/handles.Setup.PointCloud.divider*10^6*0.6;%us
    [handles.Setup,Cloud.sequenceid1D.y(i)] = function_StoreImages_DMD_timecontrol(handles.Setup, DMDFrames, illuminatetime, Stim.Npeaks);
%             [handles.Setup,Cloud.sequenceid1D.x(i)] = function_StoreImages_DMD(handles.Setup,DMDFrames);
end
for i=1:Nz
    if status == 0; disp('Procedure interrupted'); break; end;
    [DMDFrames] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string'))+UZ(i),str2num(get(handles.edit19,'string')));
    DMDFrames = 255*uint8(gather(DMDFrames));
%             DMDFrames = 255*cat(3,DMDFrames,ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY, 'uint8'));
    illuminatetime=handles.Setup.PointCloud.CycleLength/handles.Setup.PointCloud.divider*10^6*0.6;%us
    [handles.Setup,Cloud.sequenceid1D.z(i)] = function_StoreImages_DMD_timecontrol(handles.Setup, DMDFrames, illuminatetime, Stim.Npeaks);
%             [handles.Setup,Cloud.sequenceid1D.x(i)] = function_StoreImages_DMD(handles.Setup,DMDFrames);
end
[handles.Setup,Cloud.blanksequenceid] = function_StoreImages_DMD_timecontrol(handles.Setup, Stim.DMDFrame_forimage(:,:,1),...
    handles.Setup.camera.FramesPerTrigger*Stim.Exposure*10^3, 1);
% outputSingleScan(handles.Setup.Daq,Stim.Output(1,:));
% pause(20);
Images_all=zeros(Stim.ROIoffset(3),Stim.ROIoffset(4),handles.Setup.camera.FramesPerTrigger*Stim.repeatNum,Nx,'uint16');
for i=1:Nx
    if status == 0; disp('Procedure interrupted'); break; end;
    for j = 1:Stim.repeatNum %do as many repetitions as needed
        if status == 0; disp('Procedure interrupted'); break; end;
        handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.sequenceid1D.x(i));
        queueOutputData(handles.Setup.Daq,Stim.Output(1:round(Stim.Npeaks/Stim.FreqHZ*handles.Setup.Daq.Rate),:));
        startForeground(handles.Setup.Daq);

        handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.blanksequenceid);
        start(handles.Setup.camera);
        outputSingleScan(handles.Setup.Daq,[0,Stim.VoltageCa,handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1),...
            handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2), 1, 1, 0, 0]);
        currentstack=squeeze(getdata(handles.Setup.camera));
        stop(handles.Setup.camera);
        Images_all(:,:,(j-1)*handles.Setup.camera.FramesPerTrigger+1:j*handles.Setup.camera.FramesPerTrigger,i)=currentstack;
        if j==1
            Images=currentstack(:,:,2:end-2);
        else
            Images=cat(3,Images,currentstack(:,:,2:end-2));
        end
        disp(['finish ' num2str(j) '/' get(handles.edit9,'string') ' repetition']);
    end
    % plot the result--------------------------------------
    [CaTrace, dF, Mask] = function_plotCaTrace(Images, Stim);
    if i==1
        timefilter=rem(1:size(Images,3),size(Images,3)/Stim.repeatNum)<size(Images,3)/Stim.repeatNum*0.5;
        ind=find(timefilter);
    end
    score=mean(dF(ind));
    axes(handles.axes2); scatter(UX(i), score,'red','filled'); hold on; xlabel('X pixels'); ylabel('Score'); title('X PPSF');
    axes(handles.axes5);plot(CaTrace);title('flurescence of ROIs');xlabel('frames');
    axes(handles.axes3);imagesc(mean(Images,3));
    green = cat(3, zeros(size(Mask)),ones(size(Mask)),zeros(size(Mask)));hold on;h=imshow(green);hold off;set(h,'AlphaData',Mask*0.1);
    colorbar;title('image of ROIs');
    axes(handles.axes4);plot(dF);title('dF/F');xlabel('time (s)');
    DataSave.X{i} = Images;
    DataSave.SX{i} = score;
    disp(['----------finished ' num2str(i) '/' num2str(Nx) ' point']);
end
DataSave.XImages=Images_all;

Images_all=zeros(Stim.ROIoffset(3),Stim.ROIoffset(4),handles.Setup.camera.FramesPerTrigger*Stim.repeatNum,Ny,'uint16');
for i=1:Ny

    if status == 0; disp('Procedure interrupted'); break; end;
    for j = 1:Stim.repeatNum %do as many repetitions as needed
        if status == 0; disp('Procedure interrupted'); break; end;
        handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.sequenceid1D.y(i));
        queueOutputData(handles.Setup.Daq,Stim.Output(1:round(Stim.Npeaks/Stim.FreqHZ*handles.Setup.Daq.Rate),:));
        startForeground(handles.Setup.Daq);

        handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.blanksequenceid);
        start(handles.Setup.camera);
        outputSingleScan(handles.Setup.Daq,[0,Stim.VoltageCa,handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1),...
            handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2), 1, 1, 0, 0]);
        currentstack=squeeze(getdata(handles.Setup.camera));
        stop(handles.Setup.camera);
        Images_all(:,:,(j-1)*handles.Setup.camera.FramesPerTrigger+1:j*handles.Setup.camera.FramesPerTrigger,i)=currentstack;
        if j==1
            Images=currentstack(:,:,2:end-3);
        else
            Images=cat(3,Images,currentstack(:,:,2:end-3));
        end
        disp(['finish ' num2str(j) '/' get(handles.edit9,'string') ' repetition']);
    end

    % plot the result
    [CaTrace, dF, Mask] = function_plotCaTrace(Images,Stim);
    if i==1
        timefilter=rem(1:size(Images,3),size(Images,3)/Stim.repeatNum)<size(Images,3)/Stim.repeatNum*0.5;
        ind=find(timefilter);
    end
    score=mean(dF(ind));
    axes(handles.axes2); scatter(UY(i), score,'blue','filled'); hold on; xlabel('Y pixels'); ylabel('Score'); title('Y PPSF');
    axes(handles.axes5);plot(CaTrace);title('flurescence of ROIs');xlabel('frames');
    axes(handles.axes3);imagesc(mean(Images,3));
    green = cat(3, zeros(size(Mask)),ones(size(Mask)),zeros(size(Mask)));hold on;h=imshow(green);hold off;set(h,'AlphaData',Mask*0.1);
    colorbar;title('image of ROIs');
    axes(handles.axes4);plot(dF);title('dF/F');xlabel('time (s)');
    DataSave.Y{i} = Images;
    DataSave.SY{i} = score;
    disp(['----------finished ' num2str(i) '/' num2str(Ny) ' point']);
end
DataSave.YImages=Images_all;

Images_all=zeros(Stim.ROIoffset(3),Stim.ROIoffset(4),handles.Setup.camera.FramesPerTrigger*Stim.repeatNum,Nz,'uint16');
for i=1:Nz
    Images_all=zeros(Stim.ROIoffset(3),Stim.ROIoffset(4),handles.Setup.camera.FramesPerTrigger*Stim.repeatNum,Nz,'uint16');
    if status == 0; disp('Procedure interrupted'); break; end;
    for j = 1:floor(str2num(get(handles.edit9,'string'))); %do as many repetitions as needed
        if status == 0; disp('Procedure interrupted'); break; end;
        handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.sequenceid1D.z(i));
        queueOutputData(handles.Setup.Daq,Stim.Output(1:round(Stim.Npeaks/Stim.FreqHZ*handles.Setup.Daq.Rate),:));
        startForeground(handles.Setup.Daq);

        handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.blanksequenceid);
        start(handles.Setup.camera);
        outputSingleScan(handles.Setup.Daq,[0,Stim.VoltageCa,handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1),...
            handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2), 1, 1, 0, 0]);
        currentstack=squeeze(getdata(handles.Setup.camera));
        stop(handles.Setup.camera);
        Images_all(:,:,(j-1)*handles.Setup.camera.FramesPerTrigger+1:j*handles.Setup.camera.FramesPerTrigger,i)=currentstack;
        if j==1
            Images=currentstack(:,:,2:end-2);
        else
            Images=cat(3,Images,currentstack(:,:,2:end-2));
        end
        disp(['finish ' num2str(j) '/' get(handles.edit9,'string') ' repetition']);
    end

    % plot the result
    [CaTrace, dF, Mask] = function_plotCaTrace(Images,Stim);
    if i==1
        timefilter=rem(1:size(Images,3),size(Images,3)/Stim.repeatNum)<size(Images,3)/Stim.repeatNum*0.5;
        ind=find(timefilter);
    end
    score=mean(dF(ind));
    axes(handles.axes2); scatter(UZ(i), score,'green','filled'); hold on; xlabel('Z pixels'); ylabel('Score'); title('Z PPSF');
    axes(handles.axes5);plot(CaTrace);title('flurescence of ROIs');xlabel('frames');
    axes(handles.axes3);imagesc(mean(Images,3));
    green = cat(3, zeros(size(Mask)),ones(size(Mask)),zeros(size(Mask)));hold on;h=imshow(green);hold off;set(h,'AlphaData',Mask*0.1);
    colorbar;title('image of ROIs');
    axes(handles.axes4);plot(dF);title('dF/F');xlabel('time (s)');
    DataSave.Z{i} = Images;
    DataSave.SZ{i} = score;
    disp(['----------finished ' num2str(i) '/' num2str(Nz) ' point']);
end
DataSave.ZImages=Images_all;
    
handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
 
% DAQstate = handles.Setup.DAQstateZero;
% outputSingleScan(handles.Setup.Daq,DAQstate);
Stim.Camera=handles.Setup.camera;
Stim.CameraSrc=handles.Setup.src;
ToSave.type = 'Digital PPSF with camera';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.DataSave = DataSave;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave','-v7.3');
    disp('Data Saved')
end
end

% --- interleave illuminate ROIs
function pushbutton62_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = handles.Setup.DAQstateZero;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
ToSave.DMDSTATE = 'interleave stim ROIs';
ToSave.Score=[];
ToSave.Data=[];

% Generate illumination trigger----------------------------------------
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);

Stim.interleaveN = str2num(get(handles.edit67,'string')); % Number of interleave patterns
Stim.Npeaks = str2num(get(handles.edit6,'string')); % Number of blue light pulses in test
Stim.repeatNum=floor(str2num(get(handles.edit9,'string')));% number of stimulation pulses     
Stim.DurationMS = str2num(get(handles.edit8,'string'))/Stim.interleaveN; % Pulse duration in ms
Cloud.CycleLength=Stim.DurationMS/1000;
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.FreqHZ=1000/(round((1000/Stim.FreqHZ)/(Stim.DurationMS*Stim.interleaveN))*(Stim.DurationMS*Stim.interleaveN));
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.VoltageCa = str2num(get(handles.edit63,'string'));
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
handles.Setup.PointCloud.GalvoOffsetVoltage=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltageBlue;
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
if Stim.DutyCycle<0.5
    Stim.Output(:,1) =double(Stim.Subclock<0.5);
else
    Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
end
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DelayBorder = round(Stim.DelayBorder/(Cloud.CycleLength/Stim.interleaveN))*(Cloud.CycleLength/Stim.interleaveN); % to avoid start from the middle of a trigger
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
Stim.NumberCycles = floor(Stim.Totalduration/Cloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.Totalduration,Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
temp=Stim.Array;
for j = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(j-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS*Stim.interleaveN/1000+(j-1)/Stim.FreqHZ);
    for i=1:Stim.interleaveN
        temp = temp+double(Stim.UUT>Stim.DelayBorder+(i-1)*Stim.DurationMS/1000+(j-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+...
            (Stim.DurationMS-Stim.DutyCycle*Stim.DurationMS/handles.Setup.PointCloud.divider)/1000+(i-1)*Stim.DurationMS/1000+(j-1)/Stim.FreqHZ);
    end
end
Stim.Array=circshift(Stim.Array',-round(Stim.DelayBorder*0.99*handles.Setup.Daq.Rate));
Stim.Baseline = Stim.Output(:,1);
% select = Stim.UUT<handles.Setup.Scorepikes.sealtestduration;
% Stim.Output(select,8) = handles.Setup.Scorepikes.sealtestvalue;
Stim.Output(:,7)=Stim.Output(:,7).*Stim.Array;
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));

Stim.Output(:,1)=Stim.Array.*Stim.Voltage;
Stim.Output(:,8)=Stim.Array.*str2num(get(handles.edit66,'string'));

% Camera settings--------------------------------------------------------
preview='off';
Stim.Exposure = 100; %double
if strcmp(preview,'off')
    handles.Setup.src.TriggerMode = 'Trigger first'; % continue capture, Edge Trigger, Trigger first
    % handles.Setup.camera.FramesPerTrigger = Inf;
    Stim.Exposure = str2num(get(handles.edit58,'string')); %double
    triggerconfig(handles.Setup.camera, 'hardware', 'Falling edge', 'Extern');
    handles.Setup.src.Exposure = Stim.Exposure; %int32 type
    handles.Setup.src.AutoContrast = 'OFF';
    handles.Setup.src.FanSpeed = 'Low';
    handles.Setup.src.ExposeOutMode = 'All Rows'; %First Row
    handles.Setup.src.PortSpeedGain = 'Port0-Speed1-100MHz-16bit-Gain1-HDR';
    handles.Setup.src.ClearMode = 'Post-Sequence';
end

if Stim.repeatNum==1 && Stim.FreqHZ<10 % long stim
    handles.Setup.camera.FramesPerTrigger = floor((1000/Stim.FreqHZ-Stim.DurationMS)*0.99/Stim.Exposure);
elseif Stim.repeatNum>1 && Stim.FreqHZ>=10 % pulse train
    handles.Setup.camera.FramesPerTrigger = floor(Stim.DelayBorder*1.8*1000/Stim.Exposure);%600
else
    handles.Setup.camera.FramesPerTrigger = 1;
end

Stim.CropMask2=zeros(size(Stim.UUT))+double(Stim.UUT>Stim.DelayBorder).*double(Stim.UUT<Stim.DelayBorder+...
    +Stim.DurationMS*Stim.interleaveN/1000+(Stim.Npeaks-1)/Stim.FreqHZ);
Stim.Output(:,5)=circshift(ones(size(Stim.CropMask2'))-Stim.CropMask2',-round(Stim.DelayBorder*0.95*handles.Setup.Daq.Rate)); % trigger camera, 0.4
Stim.Output(1:round((Stim.DurationMS/1000+(Stim.Npeaks-1)/Stim.FreqHZ)*handles.Setup.Daq.Rate),5)=0;
Stim.Output(:,2) = Stim.VoltageCa;
Stim.Output(:,6)=circshift(temp',-round(Stim.DelayBorder*0.99*handles.Setup.Daq.Rate))+Stim.Output(:,5);

if get(handles.checkbox4,'Value') == 1
    Stim.Output(:,3)=handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1);
    Stim.Output(:,4)=handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2);
%     Stim.Output(:,3)=handles.Setup.PointCloud.GalvoOffsetVoltageBlue(1);
%     Stim.Output(:,4)=handles.Setup.PointCloud.GalvoOffsetVoltageBlue(2);
else
    Stim.galvomask=circshift(ones(size(Stim.CropMask2'))-Stim.CropMask2',-round(Stim.DelayBorder*0.99*handles.Setup.Daq.Rate));
    Stim.Output(:,3)=Stim.Output(:,3).*Stim.Array+...
        Stim.galvomask.*handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1);
    Stim.Output(:,4)=Stim.Output(:,4).*Stim.Array+...
        Stim.galvomask.*handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2);
end

[idy,idx] = find(ToSave.ImageMask);
Stim.Padarray=0;
Stim.ROIoffset=[min(idx)-Stim.Padarray min(idy)-Stim.Padarray max(idx)-min(idx)+Stim.Padarray max(idy)-min(idy)+Stim.Padarray];
handles.Setup.camera.ROIPosition = Stim.ROIoffset;

% calculate ROI pattern on DMD coordinate------------------------------------------
load('AffineMatrix.mat'); % mytform, generated by register button
DMDsize = imref2d([handles.Setup.DMD.LX,handles.Setup.DMD.LY]);
Stim.ROI = imwarp(ToSave.ImageMask,mytform,'OutputView',DMDsize);
Stim.StimROI = imwarp(ToSave.StimMask,mytform,'OutputView',DMDsize);
axes(handles.axes2); hold off; cla;axes(handles.axes5); hold off; cla;
DataSave.Images={};DataSave.Iref={};
Mask=ToSave.ImageMaskLabel(Stim.ROIoffset(2):Stim.ROIoffset(2)+Stim.ROIoffset(4)-1, ...
     Stim.ROIoffset(1):Stim.ROIoffset(1)+Stim.ROIoffset(3)-1);
Stim.ROILabel=bwlabel(Stim.StimROI, 8);
N=max(Mask(:));
for ii=1:N
    ROIarea(ii)=numel(nonzeros(Mask==ii));
end
Cloud.permutepattern=[];

if Stim.repeatNum==1 && Stim.FreqHZ<10 % long stim
    Stim.cutoffFreq=Stim.FreqHZ;
elseif Stim.repeatNum>1 && Stim.FreqHZ>=10 % pulse train
    Stim.cutoffFreq=1/(Stim.DelayBorder*2);
else
    Stim.cutoffFreq=0.5; % other condition
end

% Project patterns
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
end
function_directfeed_DMD( handles.Setup,zeros(handles.Setup.DMD.LX,handles.Setup.DMD.LY,'uint8'));
outputSingleScan(handles.Setup.Daq,Stim.Output(1,:));
disp('wait the yellow laser to be stable...');
pause(30); % wait the yellow laser to be stable;

Stim.permute='true'; 
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2301
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'slave');
end

Stim.Nregions=str2num(get(handles.edit64,'string'));% number of big regions
DMDFrameROI=zeros(handles.Setup.DMD.LX,handles.Setup.DMD.LY,Stim.Nregions*Stim.interleaveN,'uint8');

%permute
[~, Stim.permROI] = function_poissonDisc_fixgrid(ToSave.ROIcenter,80);
if size(Stim.permROI,1)~=size(ToSave.ROIcenter,1)
    disp('Poisson disk sampling error');
end
Nsubregions=floor(N/(Stim.Nregions*Stim.interleaveN));
for indcomb=1:Stim.Nregions*Stim.interleaveN
    if indcomb<Stim.Nregions*Stim.interleaveN
        indtemp=Stim.permROI(1+(indcomb-1)*Nsubregions:indcomb*Nsubregions);
    else
        indtemp=Stim.permROI(1+(indcomb-1)*Nsubregions:end);
    end
    for kk=1:length(indtemp)
        temp=Stim.ROILabel;
        temp(temp~=indtemp(kk))=0;
        temp(temp~=0)=255;
        DMDFrameROI(:,:,indcomb)=DMDFrameROI(:,:,indcomb)+uint8(temp);
    end
    figure(800);subplot(Stim.Nregions,Stim.interleaveN,indcomb);
    imagesc(DMDFrameROI(:,:,indcomb));axis image;title(num2str(indtemp));
end
       
    % compute DMD patterns of all combo
    illuminatetime=Stim.DurationMS*1000*0.9;%us
    for i=1:Stim.Nregions
        if i<Stim.Nregions
            [handles.Setup,Cloud.permutepattern(i)] = function_StoreImages_DMD_timecontrol(handles.Setup, DMDFrameROI(:,:,1+(i-1)*Stim.interleaveN:i*Stim.interleaveN), illuminatetime, Stim.Npeaks);
        else
            [handles.Setup,Cloud.permutepattern(i)] = function_StoreImages_DMD_timecontrol(handles.Setup, DMDFrameROI(:,:,1+(i-1)*Stim.interleaveN:end), illuminatetime, Stim.Npeaks);
        end
    end
    [handles.Setup,Cloud.allROIid] = function_StoreImages_DMD_timecontrol(handles.Setup, uint8(Stim.ROI)*255,...
    handles.Setup.camera.FramesPerTrigger*Stim.Exposure*1.5*10^3, 1);

    Npatterns=numel(Cloud.permutepattern);
    Stim.OutputStim=Stim.Output(1:round(Stim.Npeaks/Stim.FreqHZ*handles.Setup.Daq.Rate),:);
    Stim.OutputImage=Stim.Output(round(Stim.Npeaks/Stim.FreqHZ*handles.Setup.Daq.Rate)+1:end,:);
    
    for i = 1:Npatterns
        if status == 0; disp('Procedure interrupted'); break; end;
        for j=1:Stim.repeatNum
            if status == 0; disp('Procedure interrupted'); break; end;
            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.permutepattern(i));
            queueOutputData(handles.Setup.Daq,Stim.OutputStim);
            startForeground(handles.Setup.Daq);

            handles.Setup = function_StartProj_DMD(handles.Setup, Cloud.allROIid);
            if strcmp(preview,'off')
                start(handles.Setup.camera);
            end
            outputSingleScan(handles.Setup.Daq,[0,Stim.VoltageCa,handles.Setup.PointCloud.GalvoOffsetVoltageYellow(1),...
             handles.Setup.PointCloud.GalvoOffsetVoltageYellow(2), 1, 1, 0, 0]);
            pause(handles.Setup.camera.FramesPerTrigger*Stim.Exposure*1.5/1000);
            if strcmp(preview,'off')
                stop(handles.Setup.camera);
                currentstack=squeeze(getdata(handles.Setup.camera));
                if j==1
                    Images=currentstack(:,:,2:end-2);
                else
                    Images=cat(3,Images,currentstack(:,:,2:end-2));
                end
            end
            disp(['finish ' num2str(j) '/' get(handles.edit9,'string') ' repetition']);
        end
        if strcmp(preview,'off')
            if i==1
                Images_all=Images;
            else
                Images_all=cat(4,Images_all,Images);
            end
        end
        disp(['--------------Pattern#=' num2str(i) '/' num2str(Npatterns) ' finished!']);
    end
    
    handles.Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', handles.Setup.DMD.deviceid);
    DAQstate = handles.Setup.DAQstateZero;
    outputSingleScan(handles.Setup.Daq,DAQstate);
    %processing data------------------------------------------------------
    disp('start processing data...');
    Nframe=size(Images,3);
    Mask3D=repmat(logical(ToSave.ImageMask(Stim.ROIoffset(2):Stim.ROIoffset(2)+Stim.ROIoffset(4)-1, ...
        Stim.ROIoffset(1):Stim.ROIoffset(1)+Stim.ROIoffset(3)-1)),[1,1,Nframe,Npatterns]);
    c=nonzeros(double(Images_all).*double(Mask3D));
    Non0=reshape(c,[numel(c)/(Nframe*Npatterns),Nframe,Npatterns]);
    Ica=zeros(Nframe, Npatterns, N);
    figure(500);set(gcf,'position',[100,100,ceil(N/10)*200,100*10]);
    for ii=1:N
        if ii==1
            DataSave.Images{ii}=Non0(1:ROIarea(1),:,:);
        else
            DataSave.Images{ii}=Non0((1+sum(ROIarea(1:ii-1))):sum(ROIarea(1:ii)),:,:);
        end
        b=squeeze(mean(DataSave.Images{ii},1));
        Ica(:,:,ii)=b;
        figure(500);subplot(10,ceil(N/10),ii);plot(b(size(b,1)/2:end,:));axis tight;
        drawnow;
    end
    %compute dF/F of the brightest ROI-----------------------------
    [~,indmaxlinear]=max(Ica(:));
    [~,idy,indmax]=ind2sub(size(Ica),indmaxlinear);
    CaTrace=Ica(:,idy,indmax)-100;
    dF=zeros(size(CaTrace));
    w=round((1000/Stim.Exposure)/Stim.cutoffFreq)*0.5;
    for ii=1:numel(CaTrace)
        if ii<=w
            F0=min(CaTrace(1:ii+w));
        elseif ii>numel(CaTrace)-w
            F0=min(CaTrace(ii-w:end));
        else
            F0=min(CaTrace(ii-w:ii+w));
        end
        dF(ii)=(CaTrace(ii)-F0)/F0;
    end
    T_ds=0:Stim.Exposure:(numel(CaTrace)-1)*Stim.Exposure;
    axes(handles.axes2);plot(T_ds,CaTrace);title(['Ca2+ Intensity of Imax ROI#' num2str(indmax) ', pattern#=' num2str(idy)]);xlabel('time [ms]');ylabel('Intensity');
    axes(handles.axes5);plot(T_ds,dF);xlim([min(T_ds) max(T_ds)]);xlabel('time [ms]');ylabel('\DeltaF/F');
    title(['\DeltaF/F of of Imax ROI#' num2str(indmax) ', pattern#=' num2str(idy)]);
disp('Finish data processing!');
Stim.Camera=handles.Setup.camera;
Stim.CameraSrc=handles.Setup.src;
ToSave.DMDFrameROI=DMDFrameROI;
ToSave.DataSave=DataSave;
ToSave.DMDROI=Stim.ROI;
ToSave.ProcessedTraces=Ica;
ToSave.Stim = Stim;
ToSave.AffineMatrix=mytform;
ToSave.status = status;
handles.Setup.src.FanSpeed = 'High';
if get(handles.checkbox1,'Value') == 1
    disp('Start saving...');
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave','-v7.3');
    disp('Data Saved')
end
end


 % number of interleave patterns
function edit67_Callback(hObject, eventdata, handles)
end

% --- # of interleave patterns
function edit67_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% ---pushbutton63: project mask for preview
function pushbutton63_Callback(hObject, eventdata, handles)
global ToSave;
ToSave.DMDSTATE = 'Project ROI mask for preview';
ProjMode = function_CheckProjMode_DMD(handles.Setup);
if ProjMode==2302
    [handles.Setup]=function_StopProj_DMD(handles.Setup);
    [handles.Setup] = function_DMDProjMode(handles.Setup,'master');
end
load('AffineMatrix.mat'); % mytform, generated by register button
DMDsize = imref2d([handles.Setup.DMD.LX,handles.Setup.DMD.LY]);
Stim.ROI = imwarp(ToSave.ImageMask,mytform,'OutputView',DMDsize);
function_directfeed_DMD( handles.Setup,uint8(Stim.ROI)*255);
global DAQstate;
outputSingleScan(handles.Setup.Daq,DAQstate);
end