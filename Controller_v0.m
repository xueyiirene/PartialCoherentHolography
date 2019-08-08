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
disp('Welcome to Partially Coherenet Holography GUI v0');disp('Nidaq Starting');
handles.Setup = Function_Load_Parameters('EphyS ON');
disp('Started')
global DAQstate;
global status; status = 1; % 1 if ok to continue
global SaveID; SaveID=0;
global ToSave;
global Stage;
global Cloud;
global StageStatus; StageStatus=0;
DAQstate=[0 0 0 0 0 0];
handles.stageposition = [0 0 0];
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0); %Reset slider
set(handles.slider2,'Value',0); % Reset slider
set(handles.edit3,'string',0); %Initial stim laser Voltage
set(handles.edit4,'string',0); %Initial Imaging laser voltageset(handles.edit6,'string',int2str(10)); %Number of pulses
set(handles.edit6,'string',int2str(10)); % Number of pulses
set(handles.edit7,'string',int2str(25)); % Frequency in Hertz
set(handles.edit8,'string',int2str(3)); % Pulse duration in ms
set(handles.edit9,'string',int2str(5)); % Number of repetitions
set(handles.edit10,'string',int2str(0)); % Delay between repeitions in seconds
set(handles.edit11,'string',num2str(3.5)); %Stim laser voltage
set(handles.edit29,'string',num2str(0.3)); %Pre stimulation delay in seconds
set(handles.edit30,'string',num2str(0.85)); %Pulsed lase Duty cycle between 0 and 1
set(handles.edit12,'string',num2str(1)); %Voltage sweep start voltage
set(handles.edit13,'string',num2str(5));%Voltage sweep End voltage
set(handles.edit14,'string',num2str(10));%Voltage sweep number of steps
set(handles.edit16,'string',num2str(0)); %Targets positions pixels X
set(handles.edit17,'string',num2str(0)); %Targets positions pixels Y
set(handles.edit18,'string',num2str(-3)); %Targets positions pixels Z
set(handles.edit19,'string',num2str(8)); %Target radius
set(handles.edit20,'string',num2str(100)); % PPSFRange of measurements on X axis in pixels
set(handles.edit21,'string',num2str(0));% PPSF Range of measurements on Y axis in pixels
set(handles.edit22,'string',num2str(100));% PPSFRange of measurements on Z axis in pixels
set(handles.edit23,'string',num2str(11));% PPSF number of points on X axis
set(handles.edit24,'string',num2str(0));% PPSF number of points on Y axis
set(handles.edit25,'string',num2str(11));% PPSF number of points on Z axis
set(handles.edit26,'string','OFF');%Stage X
set(handles.edit27,'string','OFF');%Stage Y
set(handles.edit28,'string','OFF');%Stage Z
set(handles.edit39,'string',num2str(0.005));%Timesweep Min stim duration [ms]
set(handles.edit40,'string',num2str(3));%Timesweep Max stim duration [ms]
set(handles.edit41,'string',num2str(15));%Timesweep number of steps
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
set(handles.checkbox2,'Value', 1);
set(handles.checkbox6,'Value', 0);
set(handles.checkbox7,'Value', 0);
function_directfeed_DMD( handles.Setup,Stim.BlankFrame );
handles.output = hObject;
guidata(hObject, handles);
end


% --- Executes Temporal Optimization
function pushbutton21_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
DAQstate = [0 0 0 0 0 0];
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);set(handles.edit4,'string',0);
Stim.Npeaks = str2num(get(handles.edit6,'string')); % Number of blue light pulses in test
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train

if get(handles.checkbox3,'Value') == 1
    Stim.Timeramp = exp(linspace(log(str2num(get(handles.edit39,'string'))/1000),log(str2num(get(handles.edit40,'string'))/1000),floor(str2num(get(handles.edit41,'string')))));
else
    Stim.Timeramp = linspace(str2num(get(handles.edit39,'string'))/1000,str2num(get(handles.edit40,'string'))/1000,floor(str2num(get(handles.edit41,'string')))); %Voltage ramp to test how much light is needed to stim...
end
Stim.DutyCycle = str2num(get(handles.edit29,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
if handles.Setup.Scorepikes.Method == 0
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
Stim.Array = Stim.UUT';
%g = figure(2);
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
DataSave = {};
for i = 1:floor(str2num(get(handles.edit9,'string'))) % Do many repetitions
    for j = 1:numel(Stim.Timeramp)
        if status == 0; disp('Procedure interrupted');  break; end;
        Stim.Output(:,1) = 0;
        select = (Stim.Array>Stim.DelayBorder) .*(Stim.Array<Stim.DelayBorder+Stim.Timeramp(j));
        Stim.Output(select>0,1) = Stim.Voltage;
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        if get(handles.checkbox2,'Value') == 1
            handles.Setup.Scorepikes.Method = 1; 
        else
            Setup.Scorepikes.Method=0; 
        end
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
        Cell.Timescore{j,i} = score;
        DataSave{j,i} = Data;
        axes(handles.axes2); 
        scatter(Stim.Timeramp(j),Cell.Timescore{j,i},'filled','red'); hold on;xlabel('Pulse duration'); ylabel('Ephys score [AU]');title(['Calibration pass ' int2str(i) 'of' int2str(floor(str2num(get(handles.edit9,'string'))))])
        axes(handles.axes3);
        plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage'); title([ 'Score = ' num2str(floor(Cell.Timescore{j,i}*100)/100)]);
        axes(handles.axes4);
        plot(Stim.UUT,Stim.Output(:,1));xlabel('Time [s]'); ylabel('Stim laser [V]');title(['Sequence ' int2str(i) ' of ' int2str(floor(str2num(get(handles.edit9,'string'))))]);

    end
    
end
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
if get(handles.checkbox4,'Value')==1
    Cloud.AngleMagnitude=0;
    set(handles.edit43,'string',num2str(Cloud.AngleMagnitude));% galvo voltage
    if get(handles.checkbox6,'Value')==1
        handles.Setup.PointCloud.GalvoOffsetVoltage=[1.2, 0.2];
        Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;    
    elseif get(handles.checkbox7,'Value')==1
        handles.Setup.PointCloud.GalvoOffsetVoltage=[-2.1, 0.1];
        Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage; 
    else
        handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
        Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;
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
DAQstate = [0 0 0 0 0 0];
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0);
set(handles.slider2,'Value',0);
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
Stim.Npeaks = str2num(get(handles.edit6,'string')); % Number of blue light pulses in test
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox6,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[1.2, 0.2];
    Cloud.AnlgeMagnitudeOffset=[1.2,0.2];
elseif get(handles.checkbox7,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[-2.1, 0.1];
    Cloud.AnlgeMagnitudeOffset=[-2.1, 0.1];
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
    Cloud.AnlgeMagnitudeOffset=[0,0];
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
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
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox6,'Value') == 1 || get(handles.checkbox7,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
end
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
Stim.Output(end,:)=0;
for i = 1:floor(str2num(get(handles.edit9,'string')));
    if status == 0; disp('Procedure interrupted');  break; end;
    [Result.DMDFrames,Result.TotalFrame] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
    function_feed_DMD(Result.DMDFrames );
    queueOutputData(handles.Setup.Daq,Stim.Output);
    Data=startForeground(handles.Setup.Daq);
    [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
    if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
    [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
    %f = figure(1);
    axes(handles.axes2);
    plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage'); hold off;title([ 'Ephys Data , Score = ' num2str(floor(score*100)/100)]);
    axes(handles.axes3);
    plot(Stim.UUT,HPData);xlabel('Time s'); ylabel('Voltage'); title([ 'High Passed, Score = ' num2str(floor(score*100)/100)]);
    axes(handles.axes4);
    plot(Stim.UUT,Stim.Output(:,1));xlabel('Time [s]'); ylabel('Stim laser [V]');title(['Sequence ' int2str(i) ' of ' int2str(floor(str2num(get(handles.edit9,'string'))))]);
    axes(handles.axes5);
    scatter(i,score,'filled','red'), xlabel('Trial #'); ylabel('score [AU]'); hold on;
    pause(str2num(get(handles.edit10,'string')));
    DataSave{i} = Data;
    disp(['Repetition #' num2str(i) '/' get(handles.edit9,'string') ' finished!']);
end

DAQstate = [0 0 0 0 0 0];
outputSingleScan(handles.Setup.Daq,DAQstate);
ToSave.status = status;
ToSave.type = 'Undefined Stimulation cycle, rawdata';
ToSave.Stim = Stim;
ToSave.Data = DataSave;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end
end

% --- Executes Optimize voltage
function pushbutton3_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
global ToSave;
global Cloud;
global SaveID; SaveID=SaveID+1;
DAQstate = [0 0 0 0 0 0];
outputSingleScan(handles.Setup.Daq,DAQstate);
set(handles.slider3,'Value',0)
set(handles.slider2,'Value',0)
set(handles.edit3,'string',0);
set(handles.edit4,'string',0);
Stim.Npeaks = str2num(get(handles.edit6,'string')); % Number of blue light pulses in test
Stim.FreqHZ= str2num(get(handles.edit7,'string')); % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = str2num(get(handles.edit8,'string')); % Pulse duration in ms
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DutyCycle = str2num(get(handles.edit29,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox6,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[1.2, 0.2];
    Cloud.AnlgeMagnitudeOffset=[1.2,0.2];
elseif get(handles.checkbox7,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[-2.1, 0.1];
    Cloud.AnlgeMagnitudeOffset=[-2.1, 0.1];
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
    Cloud.AnlgeMagnitudeOffset=[0,0];
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
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
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox6,'Value') == 1 || get(handles.checkbox7,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
end
%g = figure(2);
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
DataSave = {};
for i = 1:floor(str2num(get(handles.edit9,'string'))) % Do many repetitions
    for j = 1:numel(Stim.Voltageramp)
        if status == 0; disp('Procedure interrupted');  break; end;
        Stim.Output(:,1) = Stim.Voltageramp(j)*Stim.Baseline.*Stim.Array';
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
        Cell.Powerscore{j,i} = score;
        DataSave{j,i} = Data;
         axes(handles.axes2); 
        scatter(Stim.Voltageramp(j),Cell.Powerscore{j,i},'filled','red'); hold on;xlabel('Voltage on laser'); ylabel('Ephys score [AU]');title(['Calibration pass ' int2str(i) 'of' int2str(floor(str2num(get(handles.edit9,'string'))))])
        axes(handles.axes3);
         plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage'); title([ 'Score = ' num2str(floor(Cell.Powerscore{j,i}*100)/100)]);
         axes(handles.axes4);
         plot(Stim.UUT,Stim.Output(:,1));xlabel('Time [s]'); ylabel('Stim laser [V]');title(['Sequence ' int2str(i) ' of ' int2str(floor(str2num(get(handles.edit9,'string'))))]);

    end
end
if status == 0
else;
    Cell.BestLaserVoltage = ginput(1); %close(g);
    scatter(Cell.BestLaserVoltage(1),Cell.BestLaserVoltage(2),'filled','blue'); title('Optimized voltage selected');
    Cell.BestLaserVoltage = Cell.BestLaserVoltage(1);
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
function_directfeed_DMD( handles.Setup,Stim.BlankFrame );
global DAQstate;
outputSingleScan(handles.Setup.Daq,DAQstate);
end

% --- Executes Disply target
function pushbutton6_Callback(hObject, eventdata, handles)
global ToSave;
global Cloud;
global DAQstate;
ToSave.DMDSTATE = 'DMD Showing one target';
handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
    if get(handles.checkbox6,'Value')==1
        handles.Setup.PointCloud.GalvoOffsetVoltage=[1.2, 0.2];
        Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;    
    elseif get(handles.checkbox7,'Value') == 1
        handles.Setup.PointCloud.GalvoOffsetVoltage=[-2.1, 0.1];
        Cloud.AnlgeMagnitudeOffset=handles.Setup.PointCloud.GalvoOffsetVoltage;
    else
        handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
        Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;
    end
DAQstate = [DAQstate(1) DAQstate(2) handles.Setup.PointCloud.GalvoOffsetVoltage(1) handles.Setup.PointCloud.GalvoOffsetVoltage(2) 0 0];
[Result.DMDFrames,Result.TotalFrame] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
function_directfeed_DMD( handles.Setup,squeeze(Result.DMDFrames(:,:,1)) );
outputSingleScan(handles.Setup.Daq,DAQstate);
end

% --- Executes   BRIGHT DMD
function pushbutton7_Callback(hObject, eventdata, handles)
global ToSave;
ToSave.DMDSTATE = 'DMD Bright';
Stim.BrightFrame = ones(handles.Setup.DMD.LX,handles.Setup.DMD.LY);
function_directfeed_DMD( handles.Setup,Stim.BrightFrame );
global DAQstate;
outputSingleScan(handles.Setup.Daq,DAQstate);
end


% --- Executes PPSF Mechanical
function pushbutton17_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
global StageStatus;
global Stage;
global Cloud;
if StageStatus==0
    Stage  = function_Start_stage( handles.Setup.MechStageComport );
    StageStatus = 1;
end
function_Zero_stage( Stage );
DAQstate = [0 0 0 0 0 0];
global ToSave;
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
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox6,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[1.2, 0.2];
    Cloud.AnlgeMagnitudeOffset=[1.2,0.2];
elseif get(handles.checkbox7,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[-2.1, 0.1];
    Cloud.AnlgeMagnitudeOffset=[-2.1, 0.1];
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
    Cloud.AnlgeMagnitudeOffset=[0,0];
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
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
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox6,'Value') == 1 || get(handles.checkbox7,'Value') == 1
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
end

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(-str2num(get(handles.edit20,'string')),str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(-str2num(get(handles.edit21,'string')),str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(-str2num(get(handles.edit22,'string')),str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;
DataSave.X = {};  DataSave.Y = {}; DataSave.Z = {}; DataSave.SX = {}; DataSave.SY = {}; DataSave.SZ = {};
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
[Result.DMDFrames,Result.TotalFrame] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
function_feed_DMD(Result.DMDFrames );
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
Stim.Output(end,:)=0;
for j = 1:floor(str2num(get(handles.edit9,'string'))); %do as many repetitions as needed
    for i = 1:floor(str2num(get(handles.edit23,'string')));
        if status == 0; disp('Procedure interrupted'); break; end;
        Position = [UX(i) 0 0];
        function_Goto_stage( Stage,Position );
        set(handles.edit26,'string',num2str(Position(1)));%Stage X
        set(handles.edit27,'string',num2str(Position(2)));%Stage Y
        set(handles.edit28,'string',num2str(Position(3)));pause(1);%Stage Z
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        [score, odata] = function_scorespikes(handles.Setup,Stim, Data );
        axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys data , Score = ' num2str(score)]);
        axes(handles.axes2);scatter(UX(i), score,'red','filled'); hold on; xlabel('X pixels'); ylabel('Score'); title('X PPSF');
        pause(str2num(get(handles.edit10,'string')));
        DataSave.X{i,j} = Data;
        DataSave.SX{i,j} = score;
    end
    for i = 1:floor(str2num(get(handles.edit24,'string')));
        if status == 0; disp('Procedure interrupted'); break; end;
        Position = [0 UY(i) 0];
        function_Goto_stage( Stage,Position );
        set(handles.edit26,'string',num2str(Position(1)));%Stage X
        set(handles.edit27,'string',num2str(Position(2)));%Stage Y
        set(handles.edit28,'string',num2str(Position(3)));pause(1);%Stage Z
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        [score, odata] = function_scorespikes(handles.Setup,Stim, Data );
        axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys data , Score = ' num2str(score)]);
        axes(handles.axes3);scatter(UY(i), score,'red','filled'); hold on;xlabel('Y pixels'); ylabel('Score'); title('Y PPSF');
        pause(str2num(get(handles.edit10,'string')));
        DataSave.Y{i,j} = Data;
        DataSave.SY{i,j} = score;
    end
    for i = 1:floor(str2num(get(handles.edit25,'string')));
        if status == 0;disp('Procedure interrupted');break; end;
        Position = [ 0 0 UZ(i)];
        function_Goto_stage( Stage,Position );
        set(handles.edit26,'string',num2str(Position(1)));%Stage X
        set(handles.edit27,'string',num2str(Position(2)));%Stage Y
        set(handles.edit28,'string',num2str(Position(3)));pause(1);%Stage Z
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
        axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys Data , Score = ' num2str(score)]);
        axes(handles.axes4); scatter(UZ(i), score,'red','filled'); hold on;xlabel('Z pixels'); ylabel('Score'); title('Z PPSF');
        pause(str2num(get(handles.edit10,'string')));
        DataSave.Z{i,j} = Data;
        DataSave.SZ{i,j} = score;
    end
end
Position = [ 0 0 0];
function_Goto_stage( Stage,Position );
set(handles.edit26,'string',num2str(Position(1)));%Stage X
set(handles.edit27,'string',num2str(Position(2)));%Stage Y
set(handles.edit28,'string',num2str(Position(3)));pause(1);%Stage Z
ToSave.type = 'Mechanical PPSF';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
    disp('Data Saved')
end
end


% --- Executes Digital PPSF with just full XZ plane
function pushbutton22_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = [0 0 0 0 0 0];
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
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox6,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[1.2, 0.2];
    Cloud.AnlgeMagnitudeOffset=[1.2,0.2];
elseif get(handles.checkbox7,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[-2.1, 0.1];
    Cloud.AnlgeMagnitudeOffset=[-2.1, 0.1];
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
    Cloud.AnlgeMagnitudeOffset=[0,0];
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup ,Cloud);
Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
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
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox6,'Value') == 1 || get(handles.checkbox7,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
end

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UZ = function_logvec( str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(-str2num(get(handles.edit20,'string')),str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UZ = linspace(-str2num(get(handles.edit22,'string')),str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UZ = UZ;
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));

[PXX,PZZ] = ndgrid(UX,UZ);
DataSave.PXX = PXX(:);
DataSave.PZZ = PZZ(:);

DataSave.XZ = {};  DataSave.SXZ = {};
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
for j = 1:floor(str2num(get(handles.edit9,'string'))) %do as many repetitions as needed
    for i = 1:numel(PXX(:))
        if status == 0; disp('Procedure interrupted'); try close(f); catch; end; break; end;
        handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
        [Result.DMDFrames,Result.TotalFrame] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string'))+DataSave.PXX(i),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string'))+DataSave.PZZ(i),str2num(get(handles.edit19,'string')));
        function_feed_DMD(Result.DMDFrames );
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
        axes(handles.axes5); plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys Data , Score = ' num2str(score)]);
        axes(handles.axes2);  scatter(DataSave.PXX(i),DataSave.PZZ(i), 25,score,'filled'); colorbar; hold on; xlabel('X pixels'); ylabel('Z pixels'); title('XZ PPSF score'); 
        pause(str2num(get(handles.edit10,'string')));
        DataSave.XZ{i,j} = Data;
        DataSave.SXZ{i,j} = score;
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
end


% --- Executes Make PPSF digitally
function pushbutton8_Callback(hObject, eventdata, handles)
global DAQstate; global status;status = 1;
DAQstate = [0 0 0 0 0 0];
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
Stim.Voltage = str2num(get(handles.edit11,'string'));
Stim.Voltageramp = linspace(str2num(get(handles.edit12,'string')),str2num(get(handles.edit13,'string')),floor(str2num(get(handles.edit14,'string')))); %Voltage ramp to test how much light is needed to stim...
Stim.DelayBorder=str2num(get(handles.edit29,'string')); % In seconds, delay before and after stim train
Stim.DutyCycle = str2num(get(handles.edit30,'string')); %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.TargetRadius = str2num(get(handles.edit19,'string')); %Radius of the target once we know location in dmd pixels
Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-handles.Setup.DMD.LX/2,handles.Setup.DMD.LX/2,handles.Setup.DMD.LX),linspace(-handles.Setup.DMD.LY/2,handles.Setup.DMD.LY/2,handles.Setup.DMD.LY));
if get(handles.checkbox6,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[1.2, 0.2];
    Cloud.AnlgeMagnitudeOffset=[1.2,0.2];
elseif get(handles.checkbox7,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[-2.1, 0.1];
    Cloud.AnlgeMagnitudeOffset=[-2.1,0.1];
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
    Cloud.AnlgeMagnitudeOffset=[0, 0];
end
[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( handles.Setup,Cloud );
Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
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
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
Stim.Baseline = Stim.Output(:,1);
if get(handles.checkbox6,'Value') == 1 || get(handles.checkbox7,'Value') == 1
    Stim.Output(:,1)=0;
    Stim.Output(:,2) = Stim.Voltage*Stim.Baseline.*Stim.Array';
else   
    Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';
    Stim.Output(:,2)=0;
end

if get(handles.checkbox3,'Value') == 1
UX = function_logvec( str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = function_logvec( str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = function_logvec( str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
else
UX = linspace(-str2num(get(handles.edit20,'string')),str2num(get(handles.edit20,'string')),floor(str2num(get(handles.edit23,'string'))));
UY = linspace(-str2num(get(handles.edit21,'string')),str2num(get(handles.edit21,'string')),floor(str2num(get(handles.edit24,'string'))));
UZ = linspace(-str2num(get(handles.edit22,'string')),str2num(get(handles.edit22,'string')),floor(str2num(get(handles.edit25,'string'))));
end
DataSave.UX = UX;
DataSave.UY = UY;
DataSave.UZ = UZ;
Cloud.CycleLength= str2num(get(handles.edit44,'string'))/1000;
Cloud.divider = floor(str2num(get(handles.edit42,'string')));
Cloud.AngleMagnitude = str2num(get(handles.edit43,'string'));
DataSave.X = {};  DataSave.Y = {}; DataSave.Z = {}; DataSave.SX = {}; DataSave.SY = {}; DataSave.SZ = {};
axes(handles.axes2); hold off; cla;axes(handles.axes3); hold off; cla;axes(handles.axes4); hold off; cla;axes(handles.axes5); hold off; cla;
Stim.Output(end,:)=0;
for j = 1:floor(str2num(get(handles.edit9,'string'))); %do as many repetitions as needed
    for i = 1:floor(str2num(get(handles.edit23,'string')));
        if status == 0; disp('Procedure interrupted'); break; end;
        handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
        [Result.DMDFrames,Result.TotalFrame] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string'))+UX(i),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
%         figure(21);imagesc(Result.TotalFrame);title(num2str(UX(i)));
        function_feed_DMD(Result.DMDFrames );
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
        axes(handles.axes5); plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys Data , Score = ' num2str(score)]);
        axes(handles.axes2);  scatter(UX(i), score,'red','filled'); hold on; xlabel('X pixels'); ylabel('Score'); title('X PPSF');
        pause(str2num(get(handles.edit10,'string')));
        DataSave.X{i,j} = Data;
        DataSave.SX{i,j} = score;
    end
    for i = 1:floor(str2num(get(handles.edit24,'string')));
        if status == 0; disp('Procedure interrupted');break; end;
        handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
        [Result.DMDFrames,Result.TotalFrame] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string'))+UY(i),str2num(get(handles.edit18,'string')),str2num(get(handles.edit19,'string')));
        function_feed_DMD(Result.DMDFrames );
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
        axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys Data , Score = ' num2str(score)]);
        axes(handles.axes3); scatter(UY(i), score,'red','filled'); hold on;xlabel('Y pixels'); ylabel('Score'); title('Y PPSF');
        pause(str2num(get(handles.edit10,'string')));
        DataSave.Y{i,j} = Data;
        DataSave.SY{i,j} = score;
    end
    for i = 1:floor(str2num(get(handles.edit25,'string')));
        if status == 0;disp('Procedure interrupted');break; end;
        handles.Setup.PointCloud.AngleMagnitude = Cloud.AngleMagnitude;
        [Result.DMDFrames,Result.TotalFrame] =  function_makespots(handles.Setup,str2num(get(handles.edit16,'string')),str2num(get(handles.edit17,'string')),str2num(get(handles.edit18,'string'))+UZ(i),str2num(get(handles.edit19,'string')));
%         figure(22);imagesc(Result.TotalFrame);title(num2str(UZ(i)));
        function_feed_DMD(Result.DMDFrames );
        queueOutputData(handles.Setup.Daq,Stim.Output);
        Data=startForeground(handles.Setup.Daq);
        [ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),handles.Setup.RCCutoffFrequencyHz );
        if get(handles.checkbox2,'Value') == 1;handles.Setup.Scorepikes.Method = 1; else Setup.Scorepikes.Method=0; end;
        [score,odata] = function_scorespikes(handles.Setup,Stim, Data );
        axes(handles.axes5);plot(Stim.UUT,odata);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys Data , Score = ' num2str(score)]);
        axes(handles.axes4); scatter(UZ(i), score,'red','filled'); hold on;xlabel('Z pixels'); ylabel('Score'); title('Z PPSF');
        pause(str2num(get(handles.edit10,'string')));
        DataSave.Z{i,j} = Data;
        DataSave.SZ{i,j} = score;        
    end
    axes(handles.axes2);cla;errorbar(UX, mean(cell2mat(DataSave.SX),2),std(cell2mat(DataSave.SX),[],2));
    axes(handles.axes3);cla;errorbar(UY, mean(cell2mat(DataSave.SY),2),std(cell2mat(DataSave.SY),[],2));
    axes(handles.axes4);cla;errorbar(UZ, mean(cell2mat(DataSave.SZ),2),std(cell2mat(DataSave.SZ),[],2));
    disp(['Repetition #' num2str(j) '/' get(handles.edit9,'string') ' finished!']);
end
try
ToSave.type = 'Digital PPSF';
ToSave.Stim = Stim;
ToSave.status = status;
ToSave.Data = DataSave;
catch;end;
if get(handles.checkbox1,'Value') == 1
    filename = [handles.Setup.SavingPath, get(handles.edit31,'string'), '_', int2str(SaveID), '_.mat'];
    save(filename,'ToSave');
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
DAQstate = [0 0 0 0 0 0];
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
DAQstate = [0 0 0 0 0 0];
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
% hObject    handle to DigitalPPSFFullXZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DAQstate; global status;status = 1;
DAQstate = [0 0 0 0 0 0];
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
DAQstate = [0 0 0 0 0 0];
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
                [ SX(i,j), ~] = function_scorespikes_averageAllpulses( handles.Setup,ToSave.Stim,cell2mat(ToSave.Data.X(i,j)));
                [ SZ(i,j), ~] = function_scorespikes_averageAllpulses( handles.Setup,ToSave.Stim,cell2mat(ToSave.Data.Z(i,j)));
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


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
global Cloud;
if get(handles.checkbox6,'Value')==1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[1.2, 0.2];
    Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;  
elseif get(handles.checkbox7,'Value') == 1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[-2.1, 0.1];
    Cloud.AnlgeMagnitudeOffset=[-2.1, 0.1];
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
    Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;
end
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Cloud;
if get(handles.checkbox6,'Value')==1
    handles.Setup.PointCloud.GalvoOffsetVoltage=[-2.1, 0.1];
    Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;    
else
    handles.Setup.PointCloud.GalvoOffsetVoltage=[0, 0];
    Cloud.AnlgeMagnitudeOffset = handles.Setup.PointCloud.GalvoOffsetVoltage;
end
end
