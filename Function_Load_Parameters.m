function [ Setup ] = Function_Load_Parameters( varargin )
%Mechanical stage parameters here :
Setup.MechStageComport = 'COM4';

% Initialize DAQ, then include Channel 1, then 2, then Galvo X, Galvo Y, DMD and Camera
% Clock
Setup.Daq = daq.createSession('ni');
Setup.Daq.Rate = 50000;
addAnalogOutputChannel(Setup.Daq,'Dev1','ao0','Voltage');
addAnalogOutputChannel(Setup.Daq,'Dev1','ao1','Voltage');
addAnalogOutputChannel(Setup.Daq,'Dev1','ao2','Voltage');
addAnalogOutputChannel(Setup.Daq,'Dev1','ao3','Voltage');
addDigitalChannel(Setup.Daq,'Dev1','Port0/Line0:1','OutputOnly');
%If Ephys is listed on, include it. 
if numel(varargin)>=1
    Setup.Ephys.Status='Ephys_ON';
    disp('The Ephys Mode is on, with Ai16 as Analog Input Channel')
    ch = addAnalogInputChannel(Setup.Daq,'Dev1','ai16','Voltage');
    ch.TerminalConfig = 'SingleEnded';
end

%Load DMD 
try
filename = 'SampleApp64bit\lib\PortabilityLayer.dll';
headname ='SampleApp64bit\include\PortabilityLayer.h';
loadlibrary(filename, headname, 'alias', 'DMD');
catch
disp('Library is already loaded')
end

Setup.SavingPath = 'C:\Users\darpa_PC\Desktop\Control_Daq\SavedData\';
Setup.DMD.LY = 1024; 
Setup.DMD.LX = 768;
Setup.Cam.LX = 1024; 
Setup.Cam.LY = 1280;

%Parameters to change to change the PSF
Setup.PointCloud.divider = 10; % How many frames pr cycle (affects galvo and DMD);10

Setup.PointCloud.phiDMD=-pi*0.24; %DMD start phase (DMD only);-0.23*pi;
Setup.spotShape='circle';

% Setup.AsymetryMatrix=[1,1;-0.0375 ,1]; % This to correct for galvo assimetrical response
Setup.AsymetryMatrix=[1,0.9;0.1 ,1];

Setup.PointCloud.phiGalvo=0; % Galvomirror phase 
Setup.PointCloud.AngleMagnitude = 2.0;
Setup.PointCloud.CycleLength = 10/1000; % in seconds
Setup.PointCloud.GalvoOffsetVoltage=[0,0];%yellow [1.2, 0.2], red [-2.1, 0.1]


%Score Spikes (only for spike recording)
Setup.RCCutoffFrequencyHz = 500; %In hertz, ephys data filter cutoff
Setup.Scorepikes.Method = 0; % set to 0 = photocurrent score, 1 = count number of spikes
Setup.Scorepikes.MinPeakDistance = 0.003;
Setup.Scorepikes.SpikeSTDDetectionthreshold = 10; %x tiems the std to be called a spike
Setup.Scorepikes.Photocurrentwindowduration = 0.010; % in seconds duration of window for scoring after onset of illumination
Setup.Scorepikes.baselineduration = 0.040; %Sampling duration for baseline
Setup.Scorepikes.sealtestduration = 0.010; %Duration of light voltagep ulse sealtest.
Setup.Scorepikes.sealtestvalue = -0.04; % in volts, sealtest 

Setup.theta=asin(1/1.33);%NA=1, water immersion OL n=1.33

% X and Z ratio from type in number to actual microns.
Setup.DataAnalysis.Xratio=1.08*0.558;%from calibration data
Setup.DataAnalysis.Zratio=1.56;
end

