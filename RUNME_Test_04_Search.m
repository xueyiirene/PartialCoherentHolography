%initialization
clear all;close all;clc
[ Setup ] = Function_Load_Parameters('EphyS ON');

Stim.lowerbound = 10; %Size of final target in pixels in the dicotomic descent in eixels
Stim.Npeaks = 10; % Number of blue light pulses in test
Stim.FreqHZ= 10; % Frequency of optogenetic stimulation in Hz
Stim.DurationMS = 10; % Pulse duration in ms
Stim.DelayBorder=0.3; % In seconds, delay before and after stim train
Stim.DutyCycle = 0.6; %Fraction 0 to 1 of illumianted Fourier circle warning for reoslution
Stim.Voltage = 1.5;
Stim.TargetRadius = 10; %Radius of the target once we know location in dmd pixels
Stim.Voltageramp = linspace(1,3,4); %Voltage ramp to test how much light is needed to stim... 


Stim.Totalduration = 2*Stim.DelayBorder + Stim.Npeaks/Stim.FreqHZ;
[Setup.XX,Setup.YY] = ndgrid(linspace(-Setup.DMD.LX/2,Setup.DMD.LX/2,Setup.DMD.LX),linspace(-Setup.DMD.LY/2,Setup.DMD.LY/2,Setup.DMD.LY));

[ Stim.UT,Stim.Output, Stim.Subclock ] = function_makeCycleClock( Setup );
Stim.Output(:,1) =double(Stim.Subclock<Stim.DutyCycle);
if Setup.Scorepikes.Method == 1
Stim.Output(:,2) = 0;
select = Stim.UT<2*Setup.Scorepikes.sealtestduration.* Stim.UT>Setup.Scorepikes.sealtestduration;
Stim.Output(select,2) = Setup.Scorepikes.sealtestvalue; 
else
Stim.Output(:,2) =0;
end
Stim.NumberCycles = floor(Stim.Totalduration/Setup.PointCloud.CycleLength);
Stim.Output = repmat(Stim.Output,[Stim.NumberCycles 1]);
[Stim.LN,Stim.LX] = size(Stim.Output);
Stim.UUT = linspace(0,Stim.NumberCycles*Stim.UT(end), Stim.LN);
Stim.Array = Stim.UUT-Stim.UUT;
for i = 1:Stim.Npeaks
    Stim.Array = Stim.Array+double(Stim.UUT>Stim.DelayBorder+(i-1)/Stim.FreqHZ).*double(Stim.UUT<Stim.DelayBorder+Stim.DurationMS/1000+(i-1)/Stim.FreqHZ);
end
%plot(Stim.Array)
Stim.Baseline = Stim.Output(:,1);
Stim.Output(:,1) = Stim.Voltage*Stim.Baseline.*Stim.Array';

%Open DMD entirely
Stim.BlankFrame = ones(Setup.DMD.LX,Setup.DMD.LY);
function_directfeed_DMD( Setup,Stim.BlankFrame );

%Yellow light on
outputSingleScan(Setup.Daq,[0 1.17 0 0 0 0]);

%Blue light on
outputSingleScan(Setup.Daq,[1.3 0 0 0 0 0]);

%Stim test
for i = 1:1
queueOutputData(Setup.Daq,Stim.Output);
Data=startForeground(Setup.Daq);
[ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),Setup.RCCutoffFrequencyHz );
score = function_scorespikes(Setup,Stim, Data );
figure(1);
subplot(3,1,1);plot(Stim.UUT,Data);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys Data , Score = ' num2str(score)]);
subplot(3,1,2);plot(Stim.UUT,HPData);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys HighPassed , Score = ' num2str(score)]);
subplot(3,1,3);plot(Stim.UUT,Stim.Output(:,1));xlabel('Time s'); ylabel('Blue laser Voltage Drive');title('Laser Stim (blue)');
pause(4);
end

%all off
outputSingleScan(Setup.Daq,[0 0 0 0 0 0]);


%Power curve
Stim.BlankFrame = ones(Setup.DMD.LX,Setup.DMD.LY);
function_directfeed_DMD( Setup,Stim.BlankFrame );
Cell.Powerscore=Stim.Voltageramp-Stim.Voltageramp;
for j = 1:numel(Stim.Voltageramp)
Stim.Output(:,1) = Stim.Voltageramp(j)*Stim.Baseline.*Stim.Array';
queueOutputData(Setup.Daq,Stim.Output);
Data=startForeground(Setup.Daq);
Cell.Powerscore(j) = function_scorespikes(Setup,Stim, Data );
end
g = figure(2);
plot(Stim.Voltageramp,Cell.Powerscore); xlabel('Voltage on laser'); ylabel('Spike count');title('Select best voltage')
Cell.BestLaserVoltage = ginput(1); close(g);
Cell.BestLaserVoltage = Cell.BestLaserVoltage(1);
disp(['The best voltage for spike testing is ' num2str(Cell.BestLaserVoltage)]);
%Update routine to the right voltage
Stim.Output(:,1) = Cell.BestLaserVoltage*Stim.Baseline.*Stim.Array';



%Descent Search
figure(1);
Cell.Window = [-Setup.DMD.LX/2,Setup.DMD.LX/2 -Setup.DMD.LY/2,Setup.DMD.LY/2];
while max(Cell.Window(2)-Cell.Window(1),Cell.Window(4)-Cell.Window(3))>Stim.lowerbound
V = function_tricut(Cell.Window);
score = [0 0 0];
Stim.displayframe = Setup.XX-Setup.XX;
for i = 1:3
frame = double(Setup.XX<V(i,2)).*double(Setup.XX>V(i,1)).*double(Setup.YY<V(i,4)).*double(Setup.YY>V(i,3));
function_directfeed_DMD( Setup,frame);
queueOutputData(Setup.Daq,Stim.Output);
Data = startForeground(Setup.Daq);
[ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),Setup.RCCutoffFrequencyHz );
score(i) = function_scorespikes(Setup,Stim, Data );
subplot(2,2,i)
plot(Stim.UUT,HPData); 
title(num2str(score(i)));
pause(0.001)
Stim.displayframe = max(Stim.displayframe,(score(i)+1)*frame);
end
subplot(2,2,4)
imagesc(Stim.displayframe);
pause(1);
[~,b] = max(score);
Cell.Window = squeeze(V(b,:));
end
%Compute_target location 
Cell.x=[(Cell.Window(1)+Cell.Window(2))/2];
Cell.y=[(Cell.Window(3)+Cell.Window(4))/2];
Cell.z =0; %Radius for defocusing
% COmpute frames and DAQ outputs for one cycle




%PPSF
Result.Positions.DX = [-200 -100 0 100 200];
Result.Positions.DY = [0 0 0 0 0];
Result.Positions.DZ = [0 0 0 0 0];
Result.Positions.PPSFSCORE = [0 0 0 0 0];
f = figure(1)
for j = 1:numel(Result.Positions.DX)
    disp(['Now computing x = ' num2str(Result.Positions.DX(j)) ', y = '  num2str(Result.Positions.DY(j)) ', z = '  num2str(Result.Positions.DZ(j))])
[Result.DMDFrames,Result.TotalFrame] =  function_makespots(Setup,Cell.x+Result.Positions.DX(j),Cell.y+Result.Positions.DY(j),Cell.z+Result.Positions.DZ(j),Stim.TargetRadius);
function_feed_DMD( Result.DMDFrames );
queueOutputData(Setup.Daq,Stim.Output);
Data=startForeground(Setup.Daq);
[ HPData, LPData ] = function_RC_filter( Data, 1/(Stim.UT(2)-Stim.UT(1)),Setup.RCCutoffFrequencyHz );
score = function_scorespikes(Setup,Stim, Data );
Result.Positions.PPSFSCORE(j) = score;
subplot(2,2,1);plot(Stim.UUT,HPData);xlabel('Time s'); ylabel('Voltage'); title([ 'Ephys Highpass Data , Score = ' num2str(score)])
subplot(2,2,3);plot(Stim.UUT,Stim.Output(:,1));xlabel('Time s'); ylabel('Blue laser Voltage Drive');title('Laser Stim (blue)')
subplot(2,2,2); imagesc(Result.TotalFrame);
subplot(2,2,4); scatter3(Result.Positions.DX(j), Result.Positions.DY(j), Result.Positions.DZ(j),50,score+1,'filled'); hold on;
end
save('PPSFData.mat','Result','Stim','Cell')
return

