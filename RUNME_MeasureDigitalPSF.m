clear all;close all;clc
%%
[ Setup ] = Function_Load_Parameters();
Test.duration = 10;%Unit: second;
STest.duration=1;
EngraveVoltage = 10; %Voltage to supply to laser for recording
VisualizeVoltage = 10;  %Voltage to supply to laser for recording
Result.Zrange = 100; % In microns, range of Zstack
Result.Xrange=50;
Result.Yrange=50;
%% Measure Digital PSF
Zstageoffset=0;
Xstageoffset=0;
Ystageoffset=0;

ZDepths = linspace(-Result.Zrange,Result.Zrange,21)+Zstageoffset;
Xranges = linspace(-Result.Xrange,Result.Xrange,11)+Xstageoffset;
Yranges = linspace(-Result.Yrange,Result.Yrange,11)+Ystageoffset;
for xx=1:numel(Xranges)
    for yy=1:numel(Yranges)
        DMDFrames_final = gpuArray(false(Setup.DMD.LX,Setup.DMD.LY,Setup.PointCloud.divider*numel(ZDepths)));
        for i=1:numel(ZDepths)
            Result.DefocusingRadius = ZDepths(i);  %Radius for defocusing
            Result.TargetRadius = 10; %Radius of the target
            DMDFrames = function_makespots(Setup,Xranges(xx),Yranges(yy),Result.DefocusingRadius,Result.TargetRadius);   
            DMDFrames_final(:,:,((i-1)*Setup.PointCloud.divider+1):(i*Setup.PointCloud.divider))=DMDFrames;
%             disp([num2str(xx) ', ' num2str(yy) ', ' num2str(i)]);
        end
        Setup= function_StoreImages_DMD(Setup, uint8(gather(DMDFrames_final))*255);
        disp(['Store current sequence #' num2str(Setup.DMD.sequenceid)]);
    end
end
%%
[ UT,Output, Subclock ] = function_makeCycleClock( Setup );
NumberCycles = 30;
 LOutput = repmat(Output,[NumberCycles 1]);
 LOutput(end,:)=0;
%%
LOutput(:,2)=0;
LOutput(1:end-1,1)=2;
Index=randperm(121)+69;
for i=1:121
    Setup.DMD.sequenceid=Index(i);
    [Setup] = function_StartProj_DMD(Setup, Setup.DMD.sequenceid,0);
    queueOutputData(Setup.Daq,LOutput);
    startForeground(Setup.Daq);
    Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', Setup.DMD.deviceid);
    if Setup.DMD.alp_returnvalue~=0
        disp('Error stop projection!');
    end
end

%%
for i=70:190
     Setup.DMD.sequenceid=i;
Setup=function_StopDMDSequence(Setup);
end

%%
function_Stop_DMD(Setup);
%%

powerblue=1.2;
powerred=0;
[ cam,Parameters ] = function_StartCam( );
ImageMax=zeros(1024/2,1280/2);
outputSingleScan(Setup.Daq,[0,0,Output(1,3),Output(1,4),0,0]);
BG=squeeze(mean(function_GetFrameCam( cam,Parameters,1),3));

for ii=1
    % Move stage
    tic;
    Sampling=size(Output,1)/10;
    Image2D=zeros(1024,1280,size(Output,1)/Sampling);%camera frame size
%     Image2D=zeros(1024,1280);
    for jj=1:size(Output,1)/Sampling
        %set voltage
        if jj==1
            outputSingleScan(Setup.Daq,[powerblue,powerred,Output(jj,3),Output(jj,4),0,0]);
            outputSingleScan(Setup.Daq,[powerblue,powerred,Output(jj,3),Output(jj,4),0,1]);
        else    
            outputSingleScan(Setup.Daq,[powerblue,powerred,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,0]);
            outputSingleScan(Setup.Daq,[powerblue,powerred,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,1]);
        end
        %get the image at this scan angle
        temp = function_GetFrameCam( cam,Parameters,1);
        Image2D(:,:,jj)=squeeze(mean(temp,3))-BG;%2D image
%         figure(2);imagesc(Image2D(:,:,jj));drawnow;caxis([0,255]);
        if jj==1
            outputSingleScan(Setup.Daq,[0,0,Output(jj,3),Output(jj,4),0,1]);
        else
            outputSingleScan(Setup.Daq,[0,0,Output((jj-1)*Sampling,3),Output((jj-1)*Sampling,4),0,1]);
        end
    end
    
    temp2D=sum(Image2D,3);% sum intensity of each plane
    ImageMax(:,:,ii)=temp2D(1024/4:1024*0.75-1,1280/4:1280*0.75-1);
    end
outputSingleScan(Setup.Daq,[0,0,0,0,0,0]);
function_StopCam(cam);


figure(1);imagesc(ImageMax);
caxis([0 max(ImageMax(:))]);
colorbar;drawnow;

% w=2;
%  [~,Indmax] = max(ImageMax(:));
%  [Indx,Indy,Indz]=ind2sub(size(ImageMax),Indmax);
%  figure(2);plot(ZDepths,squeeze(mean(mean(ImageMax(Indx-w:Indx+w,Indy-w:Indy+w,:),1),2)));
