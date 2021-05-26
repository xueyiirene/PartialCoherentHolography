clear all;close all;clc
%%
[ Setup ] = Function_Load_Parameters();
%% point scan
r=20;
x_offset=round(Setup.DMD.LX/8);
y_offset=round(Setup.DMD.LY/8);
Nx=floor((Setup.DMD.LX-2*x_offset)/r/2);
Ny=floor((Setup.DMD.LY-2*y_offset)/r/2);
DMDFrames_final = false(Setup.DMD.LX,Setup.DMD.LY,Nx*Ny);
k=1;
for xx=x_offset+(r+1:2*r:Nx*2*r)
    for yy=y_offset+(r+1:2*r:Ny*2*r)
        R = round(0.5*r*(1+abs(xx-Setup.DMD.LX/2)/(Setup.DMD.LX/2))*(1+abs(yy-Setup.DMD.LY/2)/(Setup.DMD.LY/2)));
%         R=r;
        DMDFrames_final(xx-R:xx+R,yy-R:yy+R,k)=true;
        k=k+1;
    end
end
DMDFrames_final=uint8(DMDFrames_final)*255;
% figure(1);imshow3D(DMDFrames_final);
%%
[Setup, sequenceid]= function_StoreImages_DMD(Setup, DMDFrames_final);
%% trigger signal
Setup.PointCloud.CycleLength=100/1000; %in second
Setup.Daq.Rate =5000; % sampling rate of NI DAQ
Setup.PointCloud.divider = 1; % simply project one pattern
[UT, Output, Subclock ] = function_makeCycleClock(Setup);
Output(:,3)=0; % force galvo to zero
Output(:,4)=0;
Output(:,1)=Output(:,6)*5;

NumberCycles = size(DMDFrames_final+1,3);%10 buffer triggers
LOutput = repmat(Output,[NumberCycles 1]);
LOutput(end,:)=0;
%% project patterns
ProjMode = function_CheckProjMode_DMD(Setup);
if ProjMode==2301
    Setup=function_StopProj_DMD(Setup);
    Setup= function_DMDProjMode(Setup,'slave');
end

[Setup] = function_StartProj_DMD(Setup,sequenceid);
queueOutputData(Setup.Daq,LOutput);
startForeground(Setup.Daq);


Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', Setup.DMD.deviceid);
if Setup.DMD.alp_returnvalue~=0
    disp('Error stop projection!');
end
%% brightfield
ProjMode = function_CheckProjMode_DMD(Setup);
if ProjMode==2302
    Setup=function_StopProj_DMD(Setup);
    Setup= function_DMDProjMode(Setup,'master');
end
calibration = zeros(1600,2560);
calibration(800-50:800+50,1280-50:1280+50)=1;
calibration = uint8(calibration*255);

brightfield = uint8(ones(1600, 2560)*255);
function_directfeed_DMD( Setup,brightfield);
outputSingleScan(Setup.Daq,[0 0 0 0 0 0 1]);
%% registration DMD and Camera
[x,y]=meshgrid(linspace(-400,600,11),linspace(-400,600,11));
DefocusingRadius=zeros(numel(x),1);
TargetRadius=4*abs(x(:).*y(:))/10^5+2;
FinalFrames=function_makespots_ori(Setup,x(:),y(:),DefocusingRadius,TargetRadius);
FinalFrames = uint8(gather(FinalFrames(:,:,1))).*255;
ProjMode = function_CheckProjMode_DMD(Setup);
if ProjMode==2302
    Setup=function_StopProj_DMD(Setup);
    Setup= function_DMDProjMode(Setup,'master');
end
function_directfeed_DMD( Setup,FinalFrames);
outputSingleScan(Setup.Daq,[0,0 0 0 0 0 1]);
%% line-illumination
ProjMode = function_CheckProjMode_DMD(Setup);
if ProjMode==2302
    Setup=function_StopProj_DMD(Setup);
    Setup= function_DMDProjMode(Setup,'master');
end
pattern=zeros(1600,2560);
for dx=100:100:2400
    w=99;
    pattern(800-w:800+w,dx-w:dx+w)=1;
    pattern=uint8(pattern*255);
    function_directfeed_DMD( Setup,pattern);
    outputSingleScan(Setup.Daq,[5 0 0 0 0 0 0]);
    outputSingleScan(Setup.Daq,[5 0 0 0 0 0 1]);
    pause(0.1);
end
% outputSingleScan(Setup.Daq,[0 0 0 0 0 0 0]);
%% stop DMD
outputSingleScan(Setup.Daq,[0 0 0 0 0 0 1]);
Setup=function_StopProj_DMD(Setup);
function_Stop_DMD(Setup);
%% stop sutter stage
fclose(Setup.SutterStage);