function [DMDFrames,TotalFrame] = function_makespots(Setup,x,y,DefocusingRadius,TargetRadius)
spotShape=Setup.spotShape;
voltageadjust = Setup.PointCloud.AngleMagnitude/2; %All calibratiosn were made at 2 volts and scale linearly...

DMDFrames = zeros(Setup.DMD.LX,Setup.DMD.LY,Setup.PointCloud.divider);
[XX,YY] = ndgrid(linspace(-Setup.DMD.LX/2,Setup.DMD.LX/2,Setup.DMD.LX),linspace(-Setup.DMD.LY/2,Setup.DMD.LY/2,Setup.DMD.LY));
% 
% Zunique=unique(DefocusingRadius);% return the defocus category
% for zz=1:length(Zunique)
%     Zind=find(DefocusingRadius==Zunique(zz));
%     for k=1:length(Zind)
%         for j=1:Setup.PointCloud.divider
%             DMDFrames(:,:,j)=DMDFrames(:,:,j)+double(((XX-x(Zind(k))-voltageadjust*Zunique(zz)*cos(j*2*pi/Setup.PointCloud.divider+Setup.PointCloud.phiDMD)).^2 ...
%             +(YY-y(Zind(k))-voltageadjust*Zunique(zz)*sin(j*2*pi/Setup.PointCloud.divider+Setup.PointCloud.phiDMD)).^2)<TargetRadius(Zind(k))^2); %Cdot
%             TotalFrame = max(TotalFrame ,j*DMDFrames(:,:,j));
%         end
%     end
% end
LN = numel(x);
TotalFrame = zeros(Setup.DMD.LX,Setup.DMD.LY);
switch spotShape
    case 'square'
        for k = 1:LN
            for j=1:Setup.PointCloud.divider
                CurrentFrame=zeros(Setup.DMD.LX,Setup.DMD.LY);
                deltaX=floor(voltageadjust*DefocusingRadius(k)*cos(j*2*pi/Setup.PointCloud.divider+Setup.PointCloud.phiDMD));
                deltaY=floor(voltageadjust*DefocusingRadius(k)*sin(j*2*pi/Setup.PointCloud.divider+Setup.PointCloud.phiDMD));
                CurrentFrame((Setup.DMD.LX/2+x(k)-deltaX-TargetRadius(k)):(Setup.DMD.LX/2+x(k)-deltaX+TargetRadius(k)),(Setup.DMD.LY/2+y(k)-deltaY-TargetRadius(k)):(Setup.DMD.LY/2+y(k)-deltaY+TargetRadius(k)))=1;
                DMDFrames(:,:,j)=DMDFrames(:,:,j)+CurrentFrame;
                TotalFrame = max(TotalFrame ,j*DMDFrames(:,:,j));
            end
            disp(k);
        end
    case 'circle'
        for k = 1:LN
            for j=1:Setup.PointCloud.divider
                DMDFrames(:,:,j)=DMDFrames(:,:,j)+double(((XX-x(k)-voltageadjust*DefocusingRadius(k)*cos(j*2*pi/Setup.PointCloud.divider+Setup.PointCloud.phiDMD)).^2+(YY-y(k)-voltageadjust*DefocusingRadius(k)*sin(j*2*pi/Setup.PointCloud.divider+Setup.PointCloud.phiDMD)).^2)<TargetRadius(k)^2); %Cdot
                TotalFrame = max(TotalFrame ,j*DMDFrames(:,:,j));
            end
            disp(k);
        end
end
end

