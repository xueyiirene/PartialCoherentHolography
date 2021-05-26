function [ UT,Output, Subclock ] = function_makeCycleClock( varargin )

if numel(varargin) == 1
    Setup = varargin{1};

%   Output: DAQ output voltage to 8 channels: AO0 (478 nm laser), 
%   AO1 (589nm laser), AO2 (galva), AO3 (galva), P0.0 (Prime95 camera,seal test),
%   P0.1 (DMD trigger), P0.2 AOM, AO0(on 3rd NI board, current injection) 
NT = floor(Setup.Daq.Rate*Setup.PointCloud.CycleLength);
UT = linspace(0,Setup.PointCloud.CycleLength,NT);
Subclock = mod(Setup.PointCloud.divider*UT/max(UT),1);
Subclock(end) = Subclock(end-1);
Output = zeros(NT,8); %Analog Lasr 1 and 2, Galvo X, Galvo Y, Digital Camera, DMD, AOM
% Output(:,1) = double(Subclock<0.2);
% Output(:,2) = double(Subclock<0.2);
Output(:,1) = double(Subclock<0.2);
Output(:,2) = double(Subclock<0.2);
Output(:,3:4) = Setup.PointCloud.AngleMagnitude*[cos(-2*pi*UT/Setup.PointCloud.CycleLength+Setup.PointCloud.phiGalvo)',sin(-2*pi*UT/Setup.PointCloud.CycleLength+Setup.PointCloud.phiGalvo)']*Setup.AsymetryMatrix+Setup.PointCloud.GalvoOffsetVoltage; 
Output(:,5)= double(Subclock>0.5);
Output(:,6) = double(Subclock>0.2);
Output(:,7) = double(Subclock<0.2);
else
    Setup = varargin{1};
    Cloud = varargin{2};
    NT = floor(Setup.Daq.Rate*Cloud.CycleLength);
    UT = linspace(0,Cloud.CycleLength,NT);
    Subclock = mod(Cloud.divider*UT/max(UT),1);
    % the last element in Subclock, is always 0 because its max(UT). This
    % is actually wrong for Subclock, because for some thresholding, like
    % Subclock < threshold, the last element should be not included. So we
    % manually correct Subclock last element
    Subclock(end) = Subclock(end-1);
    
    Output = zeros(NT,8); %Analog Lasr 1 and 2, Galvo X, Galvo Y, Digital Camera, DMD, AOM
    Output(:,1) = double(Subclock<0.2);
    Output(:,2) = double(Subclock<0.2);
    Output(:,3:4) = Cloud.AngleMagnitude*[cos(-2*pi*UT/Cloud.CycleLength+Setup.PointCloud.phiGalvo)',sin(-2*pi*UT/Cloud.CycleLength+Setup.PointCloud.phiGalvo)']*Setup.AsymetryMatrix+Cloud.AnlgeMagnitudeOffset; 
    Output(:,5)= double(Subclock>0.5);
    Output(:,6) = double(Subclock>0.2); 
    Output(:,7) = double(Subclock<0.2);
    disp('varargin~=1');
    
end

%To clean up and make usre all lasers are off 
% Output(end,:)=0;

end

