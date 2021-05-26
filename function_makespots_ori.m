function [DMDFrames] = function_makespots_ori(Setup,x,y,DefocusingRadius,TargetRadius)
% spotShape=Setup.spotShape;
voltageadjust = Setup.PointCloud.AngleMagnitude/2; %All calibratiosn were made at 2 volts and scale linearly...
divider=Setup.PointCloud.divider;
phi=Setup.PointCloud.phiDMD;
DMDFrames = gpuArray(zeros(Setup.DMD.LX,Setup.DMD.LY,divider));
[XX,YY] = ndgrid(linspace(-Setup.DMD.LX/2,Setup.DMD.LX/2,Setup.DMD.LX),linspace(-Setup.DMD.LY/2,Setup.DMD.LY/2,Setup.DMD.LY));
XX=gpuArray(XX);
YY=gpuArray(YY);
LN = numel(x);
% TotalFrame = gpuArray(zeros(Setup.DMD.LX,Setup.DMD.LY));
% switch spotShape
%     case 'square'
%         for k = 1:LN
%             for j=1:Setup.PointCloud.divider
%                 CurrentFrame=zeros(Setup.DMD.LX,Setup.DMD.LY);
%                 deltaX=floor(voltageadjust*DefocusingRadius(k)*cos(j*2*pi/Setup.PointCloud.divider+Setup.PointCloud.phiDMD));
%                 deltaY=floor(voltageadjust*DefocusingRadius(k)*sin(j*2*pi/Setup.PointCloud.divider+Setup.PointCloud.phiDMD));
%                 CurrentFrame((Setup.DMD.LX/2+x(k)-deltaX-TargetRadius(k)):(Setup.DMD.LX/2+x(k)-deltaX+TargetRadius(k)),(Setup.DMD.LY/2+y(k)-deltaY-TargetRadius(k)):(Setup.DMD.LY/2+y(k)-deltaY+TargetRadius(k)))=1;
%                 DMDFrames(:,:,j)=DMDFrames(:,:,j)+CurrentFrame;
%                 TotalFrame = max(TotalFrame ,j*DMDFrames(:,:,j));
%             end
%             disp(k);
%         end
%     case 'circle'
        for k = 1:LN
            for j=1:divider
                DMDFrames(:,:,j)=DMDFrames(:,:,j)+double(((XX-x(k)-voltageadjust*DefocusingRadius(k)*cos(j*2*pi/divider+phi)).^2+(YY-y(k)-voltageadjust*DefocusingRadius(k)*sin(j*2*pi/divider+phi)).^2)<TargetRadius(k)^2); %Cdot
%                 TotalFrame = max(TotalFrame ,j*DMDFrames(:,:,j));
            end
        end
% end
%DMDFrames=uint8(gather(DMDFrames))*255;
%TotalFrame=gather(TotalFrame);
% DMDFrames=logical(DMDFrames);
end

