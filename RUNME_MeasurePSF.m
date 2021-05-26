clear all;close all;clc
%%
[ Setup ] = Function_Load_Parameters();
preload=0;
Test.duration = 10;%Unit: second;
STest.duration=1;
EngraveVoltage = 10; %Voltage to supply to laser for recording
VisualizeVoltage = 10;  %Voltage to supply to laser for recording
Result.Zrange = 200; % In microns, range of Zstack
Result.Nsteps = 200; % Number of slices in Z stack
if preload==0
    Result.DefocusingRadius =0;  %Radius for defocusing
    Result.TargetRadius = 10*ones(size(Result.DefocusingRadius)); %Radius of the target
    % [x,y]=meshgrid(linspace(-300,300,5),linspace(-300,300,5));
    x=zeros(size(Result.DefocusingRadius));
    y=zeros(size(Result.DefocusingRadius));

    % Compute frames and DAQ outputs for one cycle
    % [Result.DMDFrames] =  function_makespots(Setup,x,y,Result.DefocusingRadius,Result.TargetRadius);
    [Result.DMDFrames] =  function_makespots_ori(Setup,x,y,Result.DefocusingRadius,Result.TargetRadius);
    Result.DMDFrames=uint8(gather(Result.DMDFrames))*255;
else
    pattern=load('C:\ResearchData\Compressive Sensing\2560x1600_hadamard_patterns.mat');
    Result.DMDFrames=uint8(Result.DMDFrames.hadamard_pattern_data(:,:,1))*255;
end

[ UT,Output, Subclock ] = function_makeCycleClock( Setup );

%%
Output(:,2)=Output(:,6);
NumberCycles =10;
 LOutput = repmat(Output,[NumberCycles 1]);
 LOutput(end,:)=0;
%  LOutput(:,2)=LOutput(:,2);

%% Load frames on DMD
% [ UT,Output, Subclock ] = function_makeCycleClock( Setup );
%   NumberCycles =1000;
%  LOutput = repmat(Output,[NumberCycles 1]);
%  LOutput(end,:)=0;
 
% Setup.DMD.SequenceControl.RepeatModeValue=NumberCycles;
Setup=function_feed_DMD( Setup, Result.DMDFrames);
%% Blank screen on DMD
alignflag = 0;
Setup=function_feed_DMD( Setup, 255*uint8(ones(Setup.DMD.LY,Setup.DMD.LX)));

%% initialization DMD
for i=1:4
%     outputSingleScan(Setup.Daq,[0,1.6,0.2,2.7,rem(i,2),rem(i,2),0]);
%     outputSingleScan(Setup.Daq,[1.5,0,0,0,rem(i,2),rem(i,2),0]);
    outputSingleScan(Setup.Daq,[0,4,-0.2,-1.2,rem(i,2),rem(i,2),0]);
end
%%
LOutput(:,2)=0.9;
LOutput(1:end-1,1)=0;
for i=1:NumberCycles
queueOutputData(Setup.Daq,LOutput);
startForeground(Setup.Daq);
end
%% Stop projection
[Setup]=function_StopProj_DMD(Setup);
%% Change projection mode
[Setup] = function_DMDProjMode(Setup,'master');
%%
Setup=function_StopDMDSequence(Setup);
%% Stop DMD
function_Stop_DMD(Setup);
%% Stop piezo
fclose(Setup.SutterStage);
%% Adjust Output and generate LOutput
% [ Setup ] = Function_Load_Parameters();
[ UT,Output, Subclock ] = function_makeCycleClock( Setup );
  NumberCycles =1000;
 LOutput = repmat(Output,[NumberCycles 1]);
 LOutput(end,:)=0;

%  queueOutputData(Setup.Daq,LOutput);
%  startBackground(Setup.Daq);

%% control stage: remember to TURN-OFF stage afterwards!
% [ Stage ] = function_Start_stage( Setup.MechStageComport );
% function_Zero_stage( Stage );
% %%
% position = [0 0 20]; 
% function_Goto_stage( Stage,position );

%% Read stage
% [ positionmicron ] = function_Get_stage(Stage );%thorlabs stage
xyz_um = getPosition(Setup.SutterStage);
%%
Result.Zrange=50;Nz=11;
Ystageoffset=0;
Zstageoffset=50;
position = [0 Ystageoffset Zstageoffset]; 
moveTime = moveTo(Setup.SutterStage,[xyz_um(1);xyz_um(2);xyz_um(3)+position(3)]);
% function_Goto_stage( Stage,position );%thorlabs stage

powerblue=0;
powerred=4;
[ cam,Parameters ] = function_StartCam( );
ZDepths = linspace(-Result.Zrange,Result.Zrange,Nz)+Zstageoffset;
ImageMax=zeros(1024/2,1280/2,length(ZDepths));
position = [0 0 ZDepths(1)];
moveTime = moveTo(Setup.SutterStage,[xyz_um(1)+position(1);xyz_um(2)+position(2);xyz_um(3)+position(3)]);
% function_Goto_stage( Stage,position );%thorlabs stage
pause(0.5);
outputSingleScan(Setup.Daq,[0,0,Output(1,3),Output(1,4),0,0,0]);
BG=squeeze(mean(function_GetFrameCam( cam,Parameters,1),3));
%%
Output(:,3:4) = Setup.PointCloud.AngleMagnitude*[cos(-2*pi*UT/Setup.PointCloud.CycleLength+Setup.PointCloud.phiGalvo)',sin(-2*pi*UT/Setup.PointCloud.CycleLength+Setup.PointCloud.phiGalvo)']*Setup.AsymetryMatrix+[-0.2,-1.2]; 
%%
for ii=1:length(ZDepths)
    % Move stage
    tic;
    position = [0 Ystageoffset ZDepths(ii)]; 
    moveTime = moveTo(Setup.SutterStage,[xyz_um(1)+position(1);xyz_um(2)+position(2);xyz_um(3)+position(3)]);
    pause(0.1);
%     Sampling=size(Output,1)/10;
    Sampling=1;
    Image2D=zeros(1024,1280,size(Output,1)/Sampling);%camera frame size
%     Image2D=zeros(1024,1280);
    for jj=1:size(Output,1)/Sampling
        %set voltage
%         if jj==1
%             outputSingleScan(Setup.Daq,[powerblue,powerred,Output(jj,3),Output(jj,4),0,0,0]);
            outputSingleScan(Setup.Daq,[powerblue,powerred,Output(jj,3),Output(jj,4),0,1,0]);
%         else    
% %             outputSingleScan(Setup.Daq,[powerblue,powerred,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,0,0]);
%             outputSingleScan(Setup.Daq,[powerblue,powerred,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,1,0]);
%         end
        %get the image at this scan angle
        temp = function_GetFrameCam( cam,Parameters,1);
        Image2D(:,:,jj)=squeeze(mean(temp,3))-BG;%2D image
%         figure(4);imagesc(Image2D(:,:,jj));drawnow;caxis([0,255]);
%         if jj==1
%             outputSingleScan(Setup.Daq,[0,0,Output(jj,3),Output(jj,4),0,1,0]);
%         else
%             outputSingleScan(Setup.Daq,[0,0,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,1,0]);
%         end
    end
    
    temp2D=sum(Image2D,3);% sum intensity of each plane
    ImageMax(:,:,ii)=temp2D(1024/4:1024*0.75-1,1280/4:1280*0.75-1);
    disp(['Finish: ' num2str(ii) '/' num2str(length(ZDepths))]);toc;
%     save(['scatterPSF\red_slice2\Image2D_' num2str(ii) '.mat'],'Image2D');
end
outputSingleScan(Setup.Daq,[0,0,0,0,0,0,0]);
function_StopCam(cam);
%%
for i=1:length(ZDepths)
figure();imagesc(ImageMax(:,:,i));
 title(['Z=' num2str(ZDepths(i))]);
 caxis([0 max(ImageMax(:))]);
 colorbar;drawnow;
end
%%
w=10;
 [~,Indmax] = max(ImageMax(:));
 [Indx,Indy,Indz]=ind2sub(size(ImageMax),Indmax);
 figure(2);plot(ZDepths,squeeze(mean(mean(ImageMax(Indx-w:Indx+w,Indy-w:Indy+w,:),1),2)));
 
 save('scatterPSF\ImageMax.mat','ImageMax','ZDepths');
%%
for kk=1:1000
for i=1:20:200
    outputSingleScan(Setup.Daq,[0,0.8,Output(i,3),Output(i,4),0,0,0]);
end
end
%% Prime95 exe
Zstageoffset=0;Ystageoffset=0;
[ cam,Parameters ] = function_StartCam( );
ZDepths = linspace(-Result.Zrange,Result.Zrange,11)+Zstageoffset;
ImageMax=zeros(1024,1280,length(ZDepths));
position = [0 0 ZDepths(1)];
function_Goto_stage( Stage,position );pause(0.5);
outputSingleScan(Setup.Daq,[0,0,Output(1,3),Output(1,4),0,0]);
BG=squeeze(mean(function_GetFrameCam( cam,Parameters,1),3));

for ii=1:length(ZDepths)
    % Move stage
    tic;
    position = [0 Ystageoffset ZDepths(ii)]; 
    function_Goto_stage( Stage,position );pause(0.5);
    queueOutputData(Setup.Daq,LOutput);
    startBackground(Setup.Daq);
    temp = function_GetFrameCam( cam,Parameters,1);
    Image2D=squeeze(mean(temp,3))-BG;%2D image
    figure(2);imagesc(Image2D);drawnow;caxis([0,255]);
    ImageMax(:,:,ii)= Image2D;
    disp(['Finish: ' num2str(ii) '/' num2str(length(ZDepths))]);
    toc;
end

function_StopCam(cam);

for i=1:length(ZDepths)
figure(1);imagesc(ImageMax(:,:,i));
 title(['Z=' num2str(ZDepths(i))]);
 caxis([0 max(ImageMax(:))]);
 colorbar;drawnow;
end


%%
% 
% %%
% %Here you wanna align the stage to center
% %Find substyage camera, put fluorescent slide (spray paint - ask ian) here 
% 
% % queueOutputData(Setup.Daq,LOutput);
% % startBackground(Setup.Daq);
% [ cam,Parameters ] = function_StartCam( );
% % function_previewCam( cam,Parameters ); %opens a window to show you how the hologram looks like
% % function_StopCam(cam);
% % close all;
% 
% % return;
% %Here we are centered
% 
% Output(:,2) = EngraveVoltage*double(Subclock<0.2)/10*3;
% SNumberCycles = floor(STest.duration/Setup.PointCloud.CycleLength);
% SOutput = repmat(Output,[SNumberCycles 1]);
% % SUT = linspace(0,Setup.PointCloud.CycleLength*SNumberCycles,NT*SNumberCycles);
% SOutput(end,:)=0;
% 
% ZDepths = linspace(-Result.Zrange/2,Result.Zrange/2,Result.Nsteps);
% [ Stage ] = function_Start_stage( Setup.MechStageComport );
% function_Zero_stage( Stage );
% % f = figure('Position' ,[20 20 800 800])
% % pause(0.01)
% 
% % position = [0 0 100]; function_Goto_stage( Stage,position ); pause(0.2);
% %wait(Setup.Daq);  
% %% measure 3D PSF, galvomirror only
% for i=1:450
%     if rem(i,150)==1
%         position = [0 0 ZDepths(1)]; function_Goto_stage( Stage,position );
%         temp = function_GetFrameCam( cam,Parameters,1);
%         temp1=squeeze(mean(temp,3));
%         if i==1
%             Data=temp1;
%         else
%             Data=cat(3,Data,temp1);
%         end
%     else
%         pause(0.2)
%     end
%     outputSingleScan(Setup.Daq,[0,1,LOutput(i,3),LOutput(i,4),1,1]);
%     disp(i)
% end
% %% timer
% queueOutputData(Setup.Daq,Output_P00_P01);
% startBackground(Setup.Daq);
% %%
% [ Data ] = function_GetFrameCam( cam,Parameters,5);
% wait(Setup.Daq);
% Result.Zero = squeeze(mean(mean(Data,4),3));
% imagesc(Result.Zero);
% title('Reference')
% pause(1)
% 
% % position = [0 100 0]; function_Goto_stage( Stage,position ); pause(2);
% % %wait(Setup.Daq);  
% % queueOutputData(Setup.Daq,SOutput);
% % startBackground(Setup.Daq);
% % [ Data ] = function_GetFrameCam( cam,Parameters,5);%return RGB+frames(n), 4-D stack
% % % wait(Setup.Daq);
% % Result.HundredMicrons = squeeze(mean(mean(Data,4),3));
% % imagesc(Result.HundredMicrons);
% % title('100 Microns')
% % pause(1)
% 
% %First displacement is kinda slow
% position = [0 0 ZDepths(1)]; function_Goto_stage( Stage,position ); pause(1);
% 
% for j = 1:10:Result.Nsteps %Number of Z steps
% position = [0 0 ZDepths(j)]; function_Goto_stage( Stage,position ); pause(0.2);
% %wait(Setup.Daq);  
% queueOutputData(Setup.Daq,SOutput);
% startBackground(Setup.Daq);
% [ Data ] = function_GetFrameCam(cam,Parameters,5);
% wait(Setup.Daq);
% Data = squeeze(mean(mean(Data,4),3));
% imagesc(Data);
% Result.Stack{j} = Data;
% title(int2str(j));
% drawnow
% end
% 
% Result.Zaxis = ZDepths;
% function_StopCam(cam);
% position = [0 0 0]; function_Goto_stage(Stage,position); pause(0.2);
% function_Stop_stage(Stage);
% save([TheName '.mat'],'Result','-v7.3');
% 
% [LX,LY] = size(Data);
% TheData = zeros(LX,LY,Result.Nsteps);
% for j = 1:Result.Nsteps
% TheData(:,:,j) = Result.Stack{j};
% end
% TheData = TheData-min(TheData(:));
% TheData = 255*TheData/max(TheData(:));
% 
% f = figure('Position' ,[20 20 800 800]);
% for j = 1:Result.Nsteps
% imagesc(TheData(:,:,j));
% caxis([0 255]);
% pause(0.1);
% end
% 
% saveastiff(uint8(TheData),[TheName '.tif']);
