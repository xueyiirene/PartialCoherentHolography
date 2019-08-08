function [ UT,Output, Subclock ] = function_makeCycleClock( varargin )

if numel(varargin) == 1
    Setup = varargin{1};
    
%UNTITLED2 Summary of this function goes here
%   Output: DAQ output voltage to 6 channels: AO0 (478 nm laser), 
%   AO1 (589nm laser), AO2 (galva), AO3 (galva), P0.0 (Prime95 camera),
%   P0.1 (DMD trigger)
NT = floor(Setup.Daq.Rate*Setup.PointCloud.CycleLength);
UT = linspace(0,Setup.PointCloud.CycleLength,NT);
Subclock = mod(Setup.PointCloud.divider*UT/max(UT),1);
Output = zeros(NT,6); %Analog Lasr 1 and 2, Galvo X, Galvo Y, Digital Camera, DMD
% Output(:,1) = double(Subclock<0.2);
% Output(:,2) = double(Subclock<0.2);
Output(:,1) = 2*double(Subclock<0.2);
Output(:,2) = double(Subclock<0.2);
Output(:,3:4) = Setup.PointCloud.AngleMagnitude*[cos(-2*pi*UT/Setup.PointCloud.CycleLength+Setup.PointCloud.phiGalvo)',sin(-2*pi*UT/Setup.PointCloud.CycleLength+Setup.PointCloud.phiGalvo)']*Setup.AsymetryMatrix+Setup.PointCloud.GalvoOffsetVoltage; 
Output(:,5)= double(Subclock>0.5);
Output(:,6) = double(Subclock>0.2);
else
    Setup = varargin{1};
    Cloud = varargin{2};
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
NT = floor(Setup.Daq.Rate*Cloud.CycleLength);
UT = linspace(0,Cloud.CycleLength,NT);
Subclock = mod(Cloud.divider*UT/max(UT),1);
Output = zeros(NT,6); %Analog Lasr 1 and 2, Galvo X, Galvo Y, Digital Camera, DMD
Output(:,1) = double(Subclock<0.2);
Output(:,2) = double(Subclock<0.2);
Output(:,3:4) = Cloud.AngleMagnitude*[cos(-2*pi*UT/Cloud.CycleLength+Setup.PointCloud.phiGalvo)',sin(-2*pi*UT/Cloud.CycleLength+Setup.PointCloud.phiGalvo)']*Setup.AsymetryMatrix+Cloud.AnlgeMagnitudeOffset; 
Output(:,5)= double(Subclock>0.5);
Output(:,6) = double(Subclock>0.2);    
disp('varargin~=1');
    
end

%To clean up and make usre all lasers are off 
Output(end,:)=0;

end

