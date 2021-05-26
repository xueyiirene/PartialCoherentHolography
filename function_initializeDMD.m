function [Setup] = function_initializeDMD( varargin )
%initialize DMD
try
filename = 'alp4395.dll';
headname ='alp.h';
loadlibrary(filename, headname, 'alias', 'DMD');
catch
disp('Library is already loaded')
end

Setup.DMD.devicenumber=0;
Setup.DMD.deviceid = uint32(0);
Setup.DMD.deviceidptr = libpointer('uint32Ptr', Setup.DMD.deviceid);
Setup.DMD.initflag=0;
% Allocate DMD
[Setup.DMD.alp_returnvalue, Setup.DMD.deviceid] = calllib('DMD', 'AlpDevAlloc', ...
    Setup.DMD.devicenumber, Setup.DMD.initflag, Setup.DMD.initflag);
if Setup.DMD.alp_returnvalue~=0
    disp('Error allocate DMD: DMD are turned off or already on!');
end
% trigger-in read TTL rising edge
Setup.DMD.TriggerIn=2005;
Setup.DMD.TriggerEdge=2009;%rising 2009;falling edge 2008
Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpDevControl', ...
Setup.DMD.deviceid, Setup.DMD.TriggerIn, Setup.DMD.TriggerEdge);
if Setup.DMD.alp_returnvalue~=0
   disp('Error running DMD as slave!');
   Setup.DMD.alp_returnvalue=0;
end

Setup.DMD.sequenceid = [];
Setup.DMD.TriggerMode=2300; %2300: change master(2301) or slave (2302) mode
Setup.DMD.SequenceControl.RepeatMode=2100; 
Setup.DMD.SequenceControl.RepeatModeValue=1;
Setup.DMD.SequenceControl.BitplaneMode=2103; 
Setup.DMD.SequenceControl.BitplaneModeValue=1;

Setup.DMD.LY = 2560; 
Setup.DMD.LX = 1600;
disp('DMD is initialized!');
end

