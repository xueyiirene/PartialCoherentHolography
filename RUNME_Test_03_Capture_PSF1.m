%  function_Stop_stage(Stage);

clear all;close all;clc
[ Setup ] = Function_Load_Parameters();
TheName = 'Blu_PSF_2_';

Test.duration = 10;%Unit: second;
STest.duration=1;
EngraveVoltage = 10; %Voltage to supply to laser for recording
VisualizeVoltage = 10;  %Voltage to supply to laser for recording
Result.Zrange = 200; % In microns, range of Zstack
Result.Nsteps = 200; % Number of slices in Z stack
% Zratio=50/85;%input Z Result.DefocusRadius/ZDeapths
%Here you select where to shoot by specifying x,y,z, and spot size
% Result.DefocusingRadius =linspace(-100,100,9); %0.27;  %Radius for defocusing
% Result.TargetRadius = 10*ones(1,9); %Radius of the target
% [Xind,Yind]=meshgrid(linspace(-200,200,3));
% x=Xind(:);
% y=Yind(:);
Result.DefocusingRadius = -20;  %Radius for defocusing
Result.TargetRadius = 10*ones(size(Result.DefocusingRadius)); %Radius of the target
x=zeros(size(Result.DefocusingRadius));
y=zeros(size(Result.DefocusingRadius));

% Compute frames and DAQ outputs for one cycle
[Result.DMDFrames,Result.TotalFrame] =  function_makespots(Setup,x,y,Result.DefocusingRadius,Result.TargetRadius);

%If needed to visualize
imagesc(Result.TotalFrame);

[ UT,Output, Subclock ] = function_makeCycleClock( Setup );
%Load frames on DMD
function_feed_DMD( Result.DMDFrames );

NumberCycles =500;
Output(:,2)=1;
LOutput = repmat(Output,[NumberCycles 1]);
LOutput(end,:)=0;
%% Blank Screen on DMD
BlankDMD=ones(Setup.DMD.LX,Setup.DMD.LY);
% BlankDMD(1:Setup.DMD.LX/2,:)=1;
function_feed_DMD( BlankDMD );
for i=1:4
    outputSingleScan(Setup.Daq,[0,0,0,0,rem(i,2),rem(i,2)]);
end
%% Generate the target
targetType='grid';
switch targetType
    case 'letter'
        raw=imread('C:\Users\darpa_PC\Desktop\PartialCoherent3D\Prime95versionSetup\logo\nihlogo.png');
        raw_sampling=raw(50:150,20:230);
        raw_temp=raw_sampling(1:16:end,1:16:end);
        raw_temp(raw_temp>0)=255;
        raw_temp(5,5)=0;
        TargetMask=255*ones(size(raw_temp))-double(raw_temp);
        ind=find(flipud(TargetMask));
        [x,y]=ind2sub(size(TargetMask),ind);
        x=(x-size(TargetMask,1)/2)*32;
        y=(y-size(TargetMask,2)/2)*32;
        y(1:15)=y(1:15)-25;
        y(21:end)=y(21:end)+25;
        Result.DefocusingRadius = [40*ones(1,15),zeros(1,5),-40*ones(1,12)];
        Result.TargetRadius = 15*ones(1,length(ind));
    case 'grid'
        SpotNum=1;
        [Xind,Yind]=meshgrid(linspace(0,0,sqrt(SpotNum)));
        Result.DefocusingRadius = linspace(0,0,SpotNum);
        Result.TargetRadius = 10*ones(1,SpotNum); %Radius of the target
        x=Xind(:);
        y=Yind(:);
    case 'peaks2D'
        load('coordinates.mat');
        ImageSize=1200;%Prime95
        TargetCoordinate=h;
        Result.DefocusingRadius =zeros(size(TargetCoordinate,1),1); %0.27;  %Radius for defocusing
        Result.TargetRadius = 8*ones(size(TargetCoordinate,1),1); %Radius of the target
        x=ImageSize/2-TargetCoordinate(:,1);
        y=TargetCoordinate(:,2)-ImageSize/2;
end
%% Change phiDMD and initialize
Result.DefocusingRadius =[-20];
Setup.PointCloud.phiDMD=-0.24*pi;
[Result.DMDFrames,Result.TotalFrame] =  function_makespots(Setup,x,y,Result.DefocusingRadius,Result.TargetRadius);
figure();imagesc(Result.TotalFrame);

function_feed_DMD( Result.DMDFrames );
% figure();imagesc(Result.DMDFrames(:,:,1));

% queueOutputData(Setup.Daq,LOutput);
% startBackground(Setup.Daq);
for i=1:10
    outputSingleScan(Setup.Daq,[5,0,0,0,rem(i,2),rem(i,2)]);
end

%% Adjust the AsymetryMatrix of galvomirror

Setup.AsymetryMatrix=[1,0.93;0.07,1];
[ UT,Output, Subclock ] = function_makeCycleClock( Setup );
NumberCycles =200;
Output(:,1)=1.5;
LOutput = repmat(Output,[NumberCycles 1]);
LOutput(end,:)=0;
LOutput_GalvoOnly=LOutput;
LOutput_GalvoOnly(:,5:6)=0;
queueOutputData(Setup.Daq,LOutput_GalvoOnly);
 startBackground(Setup.Daq);
%% temp galvomirror only
for i=1:10*size(Output,1)
    outputSingleScan(Setup.Daq,[1,1.5,LOutput(i,3),LOutput(i,4),0,0]);
end
%% temp DMD only
% for k=1:10
% for i=1:size(Output,1)
%     outputSingleScan(Setup.Daq,[0,1,0,0,LOutput(i,5),LOutput(i,6)]);
% end
% end
LOutput_DMDonly=LOutput;
LOutput_DMDonly(:,3:4)=0;

 queueOutputData(Setup.Daq,LOutput_DMDonly);
 startBackground(Setup.Daq);

%% temp DMD+Galvomirror, adjust DMD phi
% Sampling=size(Output,1)/10;
% for k=1:1000
% for jj=1:size(Output,1)/Sampling
%     if jj==1
%         outputSingleScan(Setup.Daq,[0,1,Output(jj,3),Output(jj,4),0,1]);
%         outputSingleScan(Setup.Daq,[0,1,Output(jj,3),Output(jj,4),0,0]);
%     else    
%         outputSingleScan(Setup.Daq,[0,1,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,1]);
%         outputSingleScan(Setup.Daq,[0,1,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,0]);
%     end
% end
% end
% Setup.AsymetryMatrix=[1,0.92;0.08,1];
% [ UT,Output, Subclock ] = function_makeCycleClock( Setup );
% NumberCycles =500;
% Output(:,2)=1;
% LOutput = repmat(Output,[NumberCycles 1]);
% LOutput(end,:)=0;

 queueOutputData(Setup.Daq,LOutput);
 startBackground(Setup.Daq);

%% Adjust Output and generate LOutput
Setup.AsymetryMatrix=[1,0.9;0.1,1];
[ Setup ] = Function_Load_Parameters();
[ UT,Output, Subclock ] = function_makeCycleClock( Setup );
 Output(1:end,5)=Output(1:end,6);
 Output(1:end-1,5)=1;
 Output(1:end,1)=Output(1:end,5)*5;
 Output(1:end,2)=Output(1:end,5)*0;
 NumberCycles =500;
 LOutput = repmat(Output,[NumberCycles 1]);
 LOutput(end,:)=0;
 
  queueOutputData(Setup.Daq,LOutput);
 startBackground(Setup.Daq);

%% piezostage only
Output(1:end,3:4)=0;
Output(end,2)=1;
lOutput = repmat(Output,[1 1]);
lOutput(end,:)=0;
 %%
Zstageoffset=0;Ystageoffset=0;
ZDepths = linspace(-Result.Zrange,Result.Zrange,101)+Zstageoffset;
position = [0 0 ZDepths(1)];
function_Goto_stage( Stage,position );
outputSingleScan(Setup.Daq,lOutput(1,:));
pause(0.5);

for ii=1:length(ZDepths)
    % Move stage
    
    tic;
    position = [0 Ystageoffset ZDepths(ii)]; 
    function_Goto_stage( Stage,position );pause(0.2);
    queueOutputData(Setup.Daq,lOutput);
    startBackground(Setup.Daq);
    toc;
    disp(['Finish: ' num2str(ii) '/' num2str(length(ZDepths))]);
end
%% control stage: remember to TURN-OFF stage afterwards!
[ Stage ] = function_Start_stage( Setup.MechStageComport );
function_Zero_stage( Stage );
% position = [0 0 0]; function_Goto_stage( Stage,position );

%% Move stage to specific Zposition
Ystageoffset=0;
Zstageoffset=0;
position = [0 Ystageoffset Zstageoffset]; function_Goto_stage( Stage,position );
%% Record defocus galvomirror scanning angle theta
[ cam,Parameters ] = function_StartCam( );
ZDepths = linspace(-Result.Zrange,Result.Zrange,20)+Zstageoffset;
ImageMax=zeros(1024,1280,length(ZDepths));
for ii=1:length(ZDepths)
    % Move stage
    position = [0 Ystageoffset ZDepths(ii)]; 
    function_Goto_stage( Stage,position );
    
    if abs(ZDepths(ii))>100
        Sampling=3;
    else
        Sampling=5;
    end
    Image2D=zeros(1024,1280,size(Output,1)/Sampling,length(ZDepths));%camera frame size
    for jj=1:size(Output,1)/Sampling
        tic;
        %set voltage
        if jj==1
            outputSingleScan(Setup.Daq,[0,1.5,Output(jj,3),Output(jj,4),0,0]);
        else    
            outputSingleScan(Setup.Daq,[0,1.5,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,0]);
        end
        %get the image at this scan angle
        temp = function_GetFrameCam( cam,Parameters,1);
        Image2D(:,:,jj,ii)=squeeze(mean(temp,3));%2D image
        toc;
    end
    temp=squeeze(Image2D(:,:,:,ii));
    ImageMax(:,:,ii)=max(temp,[],3);% max intensity of each time frame
    disp(['Finish: ' num2str(ii) '/' num2str(length(ZDepths))]);
end
function_StopCam(cam);
%% Execute this one: Record defocus spot along Z
Zstageoffset=0;Ystageoffset=0;
[ cam,Parameters ] = function_StartCam( );
ZDepths = linspace(10,50,41)+Zstageoffset;
ImageMax=zeros(1024,1280,length(ZDepths));
position = [0 0 ZDepths(1)];
function_Goto_stage( Stage,position );pause(0.5);
outputSingleScan(Setup.Daq,[0,0,Output(1,3),Output(1,4),0,0]);
BG=squeeze(mean(function_GetFrameCam( cam,Parameters,1),3));
% for ii=1:length(ZDepths)
for ii=1:length(ZDepths)
    % Move stage
    tic;
    position = [0 Ystageoffset ZDepths(ii)]; 
    function_Goto_stage( Stage,position );pause(0.5);
    Sampling=size(Output,1)/10;
    Image2D=zeros(1024,1280,size(Output,1)/Sampling,length(ZDepths));%camera frame size
    for jj=1:size(Output,1)/Sampling
        
        %set voltage
        if jj==1
            outputSingleScan(Setup.Daq,[2,0,Output(jj,3),Output(jj,4),0,0]);
            outputSingleScan(Setup.Daq,[2,0,Output(jj,3),Output(jj,4),0,1]);
        else    
            outputSingleScan(Setup.Daq,[2,0,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,0]);
            outputSingleScan(Setup.Daq,[2,0,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,1]);
        end
        %get the image at this scan angle
        temp = function_GetFrameCam( cam,Parameters,1);
        Image2D(:,:,jj,ii)=squeeze(mean(temp,3))-BG;%2D image
        
    end
    
    ImageMax(:,:,ii)=sum(Image2D(:,:,:,ii),3);% sum intensity of each plane
    disp(['Finish: ' num2str(ii) '/' num2str(length(ZDepths))]);toc;
end
function_StopCam(cam);

for i=1:length(ZDepths)
figure(1);imagesc(ImageMax(:,:,i));
 title(['Z=' num2str(ZDepths(i))]);
 caxis([0 max(ImageMax(:))]);
 colorbar;drawnow;
end


%% Prime95 exe
Zstageoffset=60;Ystageoffset=0;
ZDepths = linspace(-Result.Zrange/2,Result.Zrange/2,100)+Zstageoffset;
position = [0 0 ZDepths(1)];
function_Goto_stage( Stage,position );pause(0.5);
outputSingleScan(Setup.Daq,[0,0,Output(1,3),Output(1,4),0,0])
for ii=1:length(ZDepths)
    % Move stage
    tic;
    position = [0 Ystageoffset ZDepths(ii)]; 
    function_Goto_stage( Stage,position );pause(0.5);
    Sampling=size(Output,1)/10;
    outputSingleScan(Setup.Daq,[0,0,Output(1,3),Output(1,4),1,0]);
    for jj=1:size(Output,1)/Sampling
        %set voltage
        if jj==1
            outputSingleScan(Setup.Daq,[0,1,Output(jj,3),Output(jj,4),1,0]);
            outputSingleScan(Setup.Daq,[0,1,Output(jj,3),Output(jj,4),1,1]);
        else    
            outputSingleScan(Setup.Daq,[0,1,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),1,0]);
            outputSingleScan(Setup.Daq,[0,1,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),1,1]);
        end
    end
    outputSingleScan(Setup.Daq,[0,0,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,0]);
    disp(['Finish: ' num2str(ii) '/' num2str(length(ZDepths))]);
    toc;
end

%% Generate defocus spot
Result.targetLoc=[0,0,100];
Result.defocusXY = function_DMDapertureLocation( Setup, Result );

%%
%Here you wanna align the stage to center
%Find substyage camera, put fluorescent slide (spray paint - ask ian) here 

% queueOutputData(Setup.Daq,LOutput);
% startBackground(Setup.Daq);
[ cam,Parameters ] = function_StartCam( );
% function_previewCam( cam,Parameters ); %opens a window to show you how the hologram looks like
% function_StopCam(cam);
% close all;

% return;
%Here we are centered

Output(:,2) = EngraveVoltage*double(Subclock<0.2)/10*3;
SNumberCycles = floor(STest.duration/Setup.PointCloud.CycleLength);
SOutput = repmat(Output,[SNumberCycles 1]);
% SUT = linspace(0,Setup.PointCloud.CycleLength*SNumberCycles,NT*SNumberCycles);
SOutput(end,:)=0;

ZDepths = linspace(-Result.Zrange/2,Result.Zrange/2,Result.Nsteps);
[ Stage ] = function_Start_stage( Setup.MechStageComport );
function_Zero_stage( Stage );
% f = figure('Position' ,[20 20 800 800])
% pause(0.01)

% position = [0 0 100]; function_Goto_stage( Stage,position ); pause(0.2);
%wait(Setup.Daq);  
%% measure 3D PSF, galvomirror only
for i=1:450
    if rem(i,150)==1
        position = [0 0 ZDepths(1)]; function_Goto_stage( Stage,position );
        temp = function_GetFrameCam( cam,Parameters,1);
        temp1=squeeze(mean(temp,3));
        if i==1
            Data=temp1;
        else
            Data=cat(3,Data,temp1);
        end
    else
        pause(0.2)
    end
    outputSingleScan(Setup.Daq,[0,1,LOutput(i,3),LOutput(i,4),1,1]);
    disp(i)
end
%% timer
queueOutputData(Setup.Daq,Output_P00_P01);
startBackground(Setup.Daq);
%%
[ Data ] = function_GetFrameCam( cam,Parameters,5);
wait(Setup.Daq);
Result.Zero = squeeze(mean(mean(Data,4),3));
imagesc(Result.Zero);
title('Reference')
pause(1)

% position = [0 100 0]; function_Goto_stage( Stage,position ); pause(2);
% %wait(Setup.Daq);  
% queueOutputData(Setup.Daq,SOutput);
% startBackground(Setup.Daq);
% [ Data ] = function_GetFrameCam( cam,Parameters,5);%return RGB+frames(n), 4-D stack
% % wait(Setup.Daq);
% Result.HundredMicrons = squeeze(mean(mean(Data,4),3));
% imagesc(Result.HundredMicrons);
% title('100 Microns')
% pause(1)

%First displacement is kinda slow
position = [0 0 ZDepths(1)]; function_Goto_stage( Stage,position ); pause(1);

for j = 1:10:Result.Nsteps %Number of Z steps
position = [0 0 ZDepths(j)]; function_Goto_stage( Stage,position ); pause(0.2);
%wait(Setup.Daq);  
queueOutputData(Setup.Daq,SOutput);
startBackground(Setup.Daq);
[ Data ] = function_GetFrameCam(cam,Parameters,5);
wait(Setup.Daq);
Data = squeeze(mean(mean(Data,4),3));
imagesc(Data);
Result.Stack{j} = Data;
title(int2str(j));
drawnow
end

Result.Zaxis = ZDepths;
function_StopCam(cam);
position = [0 0 0]; function_Goto_stage(Stage,position); pause(0.2);
function_Stop_stage(Stage);
save([TheName '.mat'],'Result','-v7.3');

[LX,LY] = size(Data);
TheData = zeros(LX,LY,Result.Nsteps);
for j = 1:Result.Nsteps
TheData(:,:,j) = Result.Stack{j};
end
TheData = TheData-min(TheData(:));
TheData = 255*TheData/max(TheData(:));

f = figure('Position' ,[20 20 800 800]);
for j = 1:Result.Nsteps
imagesc(TheData(:,:,j));
caxis([0 255]);
pause(0.1);
end

saveastiff(uint8(TheData),[TheName '.tif']);
