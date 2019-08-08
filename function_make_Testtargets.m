function [ frames, Output, UT,CycleVec, ClockVec  ] = function_make_Testtargets(Setup)

%This Code generates a bunch of targets positions are specified in DMD
%pixels, centered on zero, X,Y,X< Size in a LN by 4 matrix


UTheta = linspace(0,2*pi,Setup.Angles+1);
frames = zeros(Setup.DMD.LX,Setup.DMD.LY,Setup.Angles);


UX = linspace(0,(Setup.DMD.LX-1),Setup.DMD.LX);
UY = linspace(0,(Setup.DMD.LY-1),Setup.DMD.LY);
UX = UX-mean(UX); UY = UY-mean(UY);
[XX,YY] = ndgrid(UX,UY);

for i = 1:Setup.Angles
f = figure('position', [100 100 Setup.DMD.LY Setup.DMD.LX])
axis off
text(0.0, 0.5, int2str(i), 'FontSize', 450);
pause(0.1)
frame = getframe(f);
data = frame.cdata;
data = squeeze(mean(data,3));
close(f)
frames(:,:,i)=double(imresize(data,[Setup.DMD.LX Setup.DMD.LY])<100); %Cdot
end



NT = floor(Setup.Daq.Rate*Setup.CycleLength);
UT = linspace(0,Setup.CycleLength,NT);
Output = zeros(NT,6);
CycleVec = floor((Setup.Angles)*mod(UT,Setup.CycleLength)/Setup.CycleLength)+1;
ClockVec = mod(UT,Setup.CycleLength/Setup.Angles);
disp(['Galvo : ' int2str(Setup.CycleLength*1000)  ' ms/ turn'])
disp(['Galvo : ' int2str(1/Setup.CycleLength)  ' Hz'])
disp(['DMD: ' num2str((Setup.CycleLength/(Setup.Angles))*1000)  ' ms/ turn'])
disp(['DMD : ' int2str(1/(Setup.CycleLength/(Setup.Angles)))  ' Hz'])

%Output format is Analog Blue, Orange, Galvo X, Galvo Y, Digital Camera, DMD
Output(:,1) = 0;
Output(:,2) = Setup.BaseYellowVoltage*double(ClockVec>(0.00*Setup.CycleLength/(Setup.Angles)));
Output(:,3) = Setup.GalvoVoltage*cos(UTheta(CycleVec));
Output(:,4) = Setup.GalvoVoltage*sin(UTheta(CycleVec));
Output(:,6) = ClockVec>(0.8*Setup.CycleLength/(Setup.Angles));



end

