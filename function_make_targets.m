function [ frames, Output, UT,CycleVec, ClockVec  ] = function_make_targets(Positions,Setup)

%This Code generates a bunch of targets positions are specified in DMD
%pixels, centered on zero, X,Y,X< Size in a LN by 4 matrix


UTheta = linspace(0,2*pi,Setup.Angles+1);
frames = zeros(Setup.DMD.LX,Setup.DMD.LY,Setup.Angles);
UX = linspace(0,(Setup.DMD.LX-1),Setup.DMD.LX);
UY = linspace(0,(Setup.DMD.LY-1),Setup.DMD.LY);
UX = UX-mean(UX); UY = UY-mean(UY);
[XX,YY] = ndgrid(UX,UY);

[LN,~] = size(Positions);

for j = 1:LN
for i = 1:Setup.Angles
frames(:,:,i)=max(frames(:,:,i),double(abs(sqrt((XX-Positions(j,3)*cos(UTheta(i)) - Positions(j,1)).^2+(YY-Positions(j,3)*sin(UTheta(i))- Positions(j,2)).^2))<Positions(j,4))); %Cdot
end

end


NT = floor(Setup.Daq.Rate*Setup.CycleLength);
UT = linspace(0,Setup.CycleLength,NT);
Output = zeros(NT,6);
CycleVec = linspace(0,2*pi,NT);
ClockVec = mod(UT,Setup.CycleLength/Setup.Angles);
ClockVec = [1 double(diff(ClockVec)<0)];
ClockVec(end)= 0;

disp(['Galvo : ' int2str(Setup.CycleLength*1000)  ' ms/ turn'])
disp(['Galvo : ' int2str(1/Setup.CycleLength)  ' Hz'])
disp(['DMD: ' num2str((Setup.CycleLength/(Setup.Angles))*1000)  ' ms/ turn'])
disp(['DMD : ' int2str(1/(Setup.CycleLength/(Setup.Angles)))  ' Hz'])

%Output format is Analog Blue, Orange, Galvo X, Galvo Y, Digital Camera, DMD
Output(:,1) = 0;
Output(:,2) = Setup.BaseYellowVoltage*double(ClockVec>(0.00*Setup.CycleLength/(Setup.Angles)));
Output(:,3) = Setup.GalvoVoltage*sin((CycleVec));
Output(:,4) = -Setup.GalvoVoltage*cos((CycleVec));
Output(:,6) = ClockVec;

%f = figure(1);
%plot(UT,Output(:,3));hold on ;
%plot(UT,Output(:,4));hold on ;
%plot(UT,Output(:,5));hold on ;
%scatter(UT,Output(:,6));hold on ;


end

