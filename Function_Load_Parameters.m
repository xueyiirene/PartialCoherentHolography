function [ Setup ] = Function_Load_Parameters( varargin )
%Mechanical stage parameters here :
% Setup.MechStageComport = 'COM3';

% Initialize DAQ, then include Channel 1, then 2, then Galvo X, Galvo Y, DMD and Camera
% Clock
% Setup.Daq = daq.createSession('ni');
Setup.Daq = daq.createSession('ni');
Setup.Daq.Rate = 500000;
addAnalogOutputChannel(Setup.Daq,'Dev2','ao0','Voltage');
addAnalogOutputChannel(Setup.Daq,'Dev2','ao1','Voltage');
addAnalogOutputChannel(Setup.Daq,'Dev2','ao2','Voltage');
addAnalogOutputChannel(Setup.Daq,'Dev2','ao3','Voltage');
addDigitalChannel(Setup.Daq,'Dev2','Port0/Line0:2','OutputOnly');
addAnalogOutputChannel(Setup.Daq,'Dev1','ao0','Voltage'); % for inject current
%If Ephys is listed on, include it. 
if numel(varargin)>=1
    Setup.Ephys.Status='Ephys_ON';
    disp('The Ephys Mode is on, with Ai16 as Analog Input Channel')
    ch = addAnalogInputChannel(Setup.Daq,'Dev2','ai16','Voltage');
    ch.TerminalConfig = 'SingleEnded';
end
Setup.DAQstateZero=[0,0,0,0,0,0,0,0];

%Load DMD 
try
filename = 'alp4395.dll';
headname ='alp.h';
loadlibrary(filename, headname, 'alias', 'DMD');
catch
disp('Library is already loaded')
end

Setup.DMD.devicenumber=0;
Setup.DMD.deviceid = uint32(0);
Setup.DMD.deviceidptr = libpointer('uint32Ptr', Setup.DMD.deviceid);
Setup.DMD.initflag=0;
% Allocate DMD
[Setup.DMD.alp_returnvalue, Setup.DMD.deviceid] = calllib('DMD', 'AlpDevAlloc', ...
    Setup.DMD.devicenumber, Setup.DMD.initflag, Setup.DMD.initflag);
if Setup.DMD.alp_returnvalue~=0
    disp('Error allocate DMD: DMD are turned off or already on!');
end
% trigger-in read TTL rising edge
Setup.DMD.TriggerIn=2005;
Setup.DMD.TriggerEdge=2009;%rising 2009;falling edge 2008
Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpDevControl', ...
Setup.DMD.deviceid, Setup.DMD.TriggerIn, Setup.DMD.TriggerEdge);
if Setup.DMD.alp_returnvalue~=0
   disp('Error running DMD as slave!');
   Setup.DMD.alp_returnvalue=0;
end

Setup.DMD.sequenceid = [];
Setup.DMD.TriggerMode=2300; %2300: change master(2301) or slave (2302) mode
Setup.DMD.SequenceControl.RepeatMode=2100; 
Setup.DMD.SequenceControl.RepeatModeValue=1;
Setup.DMD.SequenceControl.BitplaneMode=2103; 
Setup.DMD.SequenceControl.BitplaneModeValue=1;

Setup.SavingPath = 'C:\ResearchData\PCH V2\ControlCode\PartialCoherentHolography-master\SavedData\';
Setup.DataFolder='DataFolder\';
Setup.DMD.LY = 2560; 
Setup.DMD.LX = 1600;
Setup.Cam.LX = 1024; 
Setup.Cam.LY = 1280;

Setup.SutterStage=sutterMP285('COM5');
%Parameters to change to change the PSF
Setup.PointCloud.divider = 10; % How many frames pr cycle (affects galvo and DMD);10

Setup.PointCloud.phiDMD=0; %DMD start phase (DMD only);-0.23*pi;
Setup.spotShape='circle';

% Setup.AsymetryMatrix=[1,1;-0.0375 ,1]; % This to correct for galvo assimetrical response
% Setup.AsymetryMatrix=[1,0.8;0,1.2];
Setup.AsymetryMatrix=[1,1;0,1.2];

Setup.PointCloud.phiGalvo=-0.49*pi; % Galvomirror phase0.5
Setup.PointCloud.AngleMagnitude = 2.0;
Setup.PointCloud.CycleLength = 4/1000; % in seconds
Setup.PointCloud.GalvoOffsetVoltageBlue=[0.8,0.5];
Setup.PointCloud.GalvoOffsetVoltageYellow = [1, 3];
Setup.PointCloud.GalvoOffsetVoltageRed=[0.6, -0.6];
Setup.PointCloud.GalvoOffsetVoltage=Setup.PointCloud.GalvoOffsetVoltageBlue;

% initialize Prime 95B camera
try
    Setup.camera = videoinput('pmimaq_2017b', 1, 'PM-Cam 1200x1200');
    Setup.src = getselectedsource(Setup.camera);
    Setup.src.ClearMode = 'Post-Sequence';
%     Setup.src.PortSpeedGain = 'Port0-Speed0-200MHz-12bit-Gain1-Full well';
    Setup.src.PortSpeedGain = 'Port0-Speed1-100MHz-16bit-Gain1-HDR';
    Setup.src.TriggerMode = 'Internal Trigger';
    % Setup.src.TriggerMode = 'Edge Trigger';
    Setup.src.Exposure = 10;
end

%Score Spikes (only for spike recording)
Setup.RCCutoffFrequencyHz = 500; %In hertz, ephys data filter cutoff
Setup.Scorepikes.Method = 0; % set to 0 = photocurrent score, 1 = count number of spikes
Setup.Scorepikes.MinPeakDistance = 0.003;
Setup.Scorepikes.SpikeSTDDetectionthreshold = 7; %x tiems the std to be called a spike
Setup.Scorepikes.Photocurrentwindowduration = 0.010; % in seconds duration of window for scoring after onset of illumination
Setup.Scorepikes.baselineduration = 0.0001; %Sampling duration for baseline,second
Setup.Scorepikes.sealtestduration = 0.050; %Duration of light voltage pulse sealtest.
Setup.Scorepikes.sealtestvalue = 0.005; % in volts, sealtest analog output

Setup.theta=asin(1/1.33);%NA=1, water immersion OL n=1.33

% X and Z ratio from type in number to actual microns.
Setup.DataAnalysis.Xratio=0.676;%from calibration data
Setup.DataAnalysis.Zratio=1.223;%0.888
end

