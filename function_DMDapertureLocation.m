function [defocusXY ] = function_DMDapertureLocation( Setup, Result )
%function_DMDapertureLocation calculates the aperture canter location on
%DMD corresponding to the defocused input location
%   Detailed explanation goes here
Z0=Result.targetLoc(3);
NT=floor(Setup.Daq.Rate*Setup.PointCloud.CycleLength);
phi=tan(linspace(0+Setup.PointCloud.phiGalvo,2*pi+Setup.PointCloud.phiGalvo,NT))';
sol_dx=sqrt((Z0*tan(Setup.theta))^2*ones(size(phi))./(1+phi.^2));
sol_dy=phi.*sol_dx;
defocusXY=[Result.targetLoc(1)+sol_dx,Result.targetLoc(2)+sol_dy];
end

