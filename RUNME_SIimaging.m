clear all;close all;clc
%%
Setup = Function_Load_Parameters();
% Setup = function_StopProj_DMD(Setup);
% Setup = function_DMDProjMode(Setup,'slave');
%% generate calibration matrix
[x,y]=meshgrid(linspace(-400,600,11),linspace(-400,600,11));
DefocusingRadius=zeros(numel(x),1);
TargetRadius=4*abs(x(:).*y(:))/10^5+2;
FinalFrames=function_makespots_ori(Setup,x(:),y(:),DefocusingRadius,TargetRadius);
%%
load('C:\ResearchData\Compressive Sensing\DMDpattern_r20_nPts50.mat');
FinalFrames=randomPattern;
%%
[ UT,Output, Subclock ] = function_makeCycleClock( Setup );
NumberCycles =50;
LOutput = repmat(Output,[NumberCycles 1]);
LOutput(:,1)=1;
LOutput(:,3:4)=0;%galvo=0
LOutput(end,:)=0;

[Setup, sequenceid]= function_StoreImages_DMD(Setup, uint8(FinalFrames(:,:,1))*255);
%%
[Setup] = function_StartProj_DMD(Setup,sequenceid);
queueOutputData(Setup.Daq,LOutput);
startForeground(Setup.Daq);
%% widefield to check sample
% Setup = function_StopProj_DMD(Setup);
% Setup = function_DMDProjMode(Setup,'master');
% FinalFrames=function_makespots_ori(Setup,0,0,0,800);
Setup=function_feed_DMD( Setup,FinalFrames);
outputSingleScan(Setup.Daq,[0,0,0,0,0,0,0]);
%%
for i=1:5
     outputSingleScan(Setup.Daq,[1,0,0,0,0,1,0]);
     pause(1);
     outputSingleScan(Setup.Daq,[1,0,0,0,0,0,0]);
     pause(1);
end

%%
% figure();imagesc(img_000000000_Default_000);
threshold=65500;
[FinalFrames] = function_generateTargetIllu(img_000000000_Default_000,Mask,t1,threshold);
figure();imagesc(FinalFrames);
%% Z scan--doesn't work T_T
Setup = function_StopProj_DMD(Setup);
Setup = function_DMDProjMode(Setup,'slave');
blank=uint8(ones(1600,2560))*255;
% FinalFrames=repmat(blank,[1,1,50]);
CurrentXYZ = getPosition(Setup.SutterStage);
dz=500;
for z=1:50
    dz=-160;
    moveTime = moveTo(Setup.SutterStage,[CurrentXYZ(1);CurrentXYZ(2);CurrentXYZ(3)+dz]);
    pause(0.1);
    Setup=function_feed_DMD( Setup,blank);
    outputSingleScan(Setup.Daq,[1,0,0,0,0,1,0]);
    outputSingleScan(Setup.Daq,[0,0,0,0,0,0,0]);
    pause(0.1);
end
%% Stop and close
outputSingleScan(Setup.Daq,[0,0,0,0,0,0,0]);
[Setup]=function_StopProj_DMD(Setup);
function_Stop_DMD(Setup);
fclose(Setup.SutterStage);
