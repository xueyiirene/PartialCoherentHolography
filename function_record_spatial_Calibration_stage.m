function [ DMDPixelPerMicron,DMDPixelperCamerapixel ] =function_record_spatial_Calibration_stage( Setup, Stage )

radius = 15;
frames = zeros(Setup.DMD.LX,Setup.DMD.LY,2);
UX = linspace(1,Setup.DMD.LX,Setup.DMD.LX);
UY = linspace(1,Setup.DMD.LY,Setup.DMD.LY);
UX = UX-mean(UX); UY = UY-mean(UY);
[XX,YY] = ndgrid(UX,UY);
frames(:,:,1)=double(abs(sqrt((XX-100).^2+YY.^2))<radius); %Cdot
frames(:,:,2)=double(abs(sqrt(XX.^2+YY.^2))<radius); %Cdot
function_feed_DMD( frames );
%Lock in this frame
Output = zeros(1,6);
outputSingleScan(Setup.Daq,Output);
Output(1,6) = 1;
outputSingleScan(Setup.Daq,Output);
Output(1,6) = 0;
outputSingleScan(Setup.Daq,Output);
[ cam ,parameters] = function_StartCam( );

%function_previewCam( cam,parameters );

Output(:,2) = Setup.YellowVoltageRef;
outputSingleScan(Setup.Daq,Output);
[ Data ] = function_GetFrameCam( cam,parameters,5);
Data = squeeze(mean(double(Data),4)); Reference.Cent = squeeze(Data(:,:,1));
Output = zeros(1,6);
outputSingleScan(Setup.Daq,Output);
Output(1,6) = 1;
outputSingleScan(Setup.Daq,Output);
Output(1,6) = 0;
outputSingleScan(Setup.Daq,Output);
Output(:,2) = Setup.YellowVoltageRef;
outputSingleScan(Setup.Daq,Output);


%Get a 50 microns frame
position = [0 0 0]; function_Goto_stage( Stage,position ); pause(0.2);
[ Data ] = function_GetFrameCam( cam,parameters,5);
Data = squeeze(mean(double(Data),4)); Reference.Zero = squeeze(Data(:,:,1));
%imagesc(Data(:,:,1)); caxis([0 255]) ; colormap gray ; title([ 'Frame zero']);
position = [0 50 0]; function_Goto_stage( Stage,position ); pause(0.2);
[ Data ] = function_GetFrameCam( cam,parameters,5);
Data = squeeze(mean(double(Data),4)); Reference.onemm = squeeze(Data(:,:,1));
%imagesc(Data(:,:,1)); caxis([0 255]) ; colormap gray ; title([ 'Frame zero']);
position = [0 0 0]; function_Goto_stage( Stage,position ); pause(0.2);
Output(1,:) = 0;
outputSingleScan(Setup.Daq,Output);
f = figure(1);
[a,b] = max(max(Reference.Zero));
[a,c] = max(Reference.Zero(:,b));
Reference.Zero(c,b);
[a,bb] = max(max(Reference.onemm));
[a,cc] = max(Reference.onemm(:,bb));
Reference.Zero(cc,bb);
pixdist = sqrt((b-bb)^2+(c-cc)^2);
Reference.pixelpermicron = pixdist/50;
subplot(2,2,1);
imagesc(Reference.Zero); hold on; scatter(b,c,'red'); axis image;
subplot(2,2,2);
imagesc(Reference.onemm); hold on; scatter(bb,cc,'red'); axis image;
[a,b] = max(max(Reference.Zero));
[a,c] = max(Reference.Zero(:,b));
Reference.Zero(c,b);
[a,bb] = max(max(Reference.Cent));
[a,cc] = max(Reference.onemm(:,bb));
Reference.Zero(cc,bb);
pixdist = sqrt((b-bb)^2+(c-cc)^2);
Reference.DMDpixelperCamPixel = pixdist/100;
subplot(2,2,3);
imagesc(Reference.Zero); hold on; scatter(b,c,'red'); axis image;
subplot(2,2,4);
imagesc(Reference.onemm); hold on; scatter(bb,cc,'red'); axis image;
pause(2);


close(f);
DMDPixelperCamerapixel = Reference.DMDpixelperCamPixel;
DMDPixelPerMicron = Reference.pixelpermicron;
function_StopCam( cam );
end

