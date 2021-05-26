function [DMDFrames] = function_makespots(Setup,x,y,DefocusingRadius,TargetRadius)
voltageadjust = Setup.PointCloud.AngleMagnitude/2; %All calibratiosn were made at 2 volts and scale linearly...
divider=Setup.PointCloud.divider;
phi=Setup.PointCloud.phiDMD;
[XX,YY] = ndgrid(gpuArray(linspace(-Setup.DMD.LX/2,Setup.DMD.LX/2,Setup.DMD.LX)),gpuArray(linspace(-Setup.DMD.LY/2,Setup.DMD.LY/2,Setup.DMD.LY)));
Ind=1:divider;
CenterX=x+voltageadjust*DefocusingRadius*cos(Ind*2*pi/divider+phi);
CenterY=y+voltageadjust*DefocusingRadius*sin(Ind*2*pi/divider+phi);
DMDFrames = gpuArray(false(Setup.DMD.LX,Setup.DMD.LY,divider));
for j=1:divider
    DMDFrames(:,:,j)=((XX-CenterX(j)).^2+(YY-CenterY(j)).^2)<TargetRadius^2; %Cdot
end
end

