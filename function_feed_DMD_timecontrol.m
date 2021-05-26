function [Setup] = function_feed_DMD_timecontrol(Setup, sequenceid, stack, illuminatetime)
% stack: 2560x1600xframes binary data
Setup.DMD.alp_returnvalue=0;
if size(stack,1)~=Setup.DMD.LY || size(stack,2)~=Setup.DMD.LX
    stack = permute(stack, [2 1 3]);
end
%allocate sequences
bitdepth=1;
picnum=size(stack,3);
sequenceidptr = libpointer('uint32Ptr', sequenceid);
[Setup.DMD.alp_returnvalue, sequenceid] = calllib('DMD', 'AlpSeqAlloc', ...
    Setup.DMD.deviceid, bitdepth, picnum, sequenceidptr);
if Setup.DMD.alp_returnvalue~=0
    disp('Error allocate sequence!');
    Setup.DMD.alp_returnvalue=0;
end
if isempty(Setup.DMD.sequenceid)
    Setup.DMD.sequenceid=sequenceid;
else
    Setup.DMD.sequenceid=cat(1,Setup.DMD.sequenceid,sequenceid);
end
%inquire sequences
% uservar = int32(0);
% uservarptr = libpointer('int32Ptr', uservar);
% Setup.DMD.inquiretype=int32(2110);
% [Setup.DMD.alp_returnvalue, uservar] = calllib('DMD', 'AlpSeqInquire', ...
%     Setup.DMD.deviceid, Setup.DMD.sequenceid, Setup.DMD.inquiretype, uservarptr);
% if Setup.DMD.alp_returnvalue~=0
%     disp('Error inquire sequence!');
%     Setup.DMD.alp_returnvalue=0;
% end
%% Stop last projection
% Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', Setup.DMD.deviceid);
%     if Setup.DMD.alp_returnvalue~=0
%         disp('Error stop projection!');
%         Setup.DMD.alp_returnvalue=0;
%     end
    
%% load frames
userarrayptr = libpointer('voidPtr', stack);%userarray: image to upload to DMD, 2560x1600 unit8 type 2D matrix
picoffset=0;
picload=size(stack,3);
Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpSeqPut', Setup.DMD.deviceid, ...
    sequenceid, picoffset, picload, userarrayptr);
if Setup.DMD.alp_returnvalue~=0
    disp('Error load sequence!');
    Setup.DMD.alp_returnvalue=0;
else
    disp([num2str(picload) ' frames loaded!']);
end

%% Repeat sequence
Setup.DMD.alp_returevalue = calllib('DMD','AlpSeqControl',Setup.DMD.deviceid,...
    sequenceid, Setup.DMD.SequenceControl.RepeatMode,...
    Setup.DMD.SequenceControl.RepeatModeValue);
if Setup.DMD.alp_returnvalue~=0
    disp('Error repeat sequence!');
    Setup.DMD.alp_returnvalue=0;
end
%% bitplane
Setup.DMD.alp_returevalue = calllib('DMD','AlpSeqControl',Setup.DMD.deviceid,...
   sequenceid, Setup.DMD.SequenceControl.BitplaneMode,...
    Setup.DMD.SequenceControl.BitplaneModeValue);
if Setup.DMD.alp_returnvalue~=0
    disp('Error setting bitplane!');
    Setup.DMD.alp_returnvalue=0;
end
%% sequence timing
Setup.DMD.illuminatetime=illuminatetime; %us Setup.PointCloud.CycleLength/Setup.PointCloud.divider*10^6*0.6;%us
picturetime=0;
triggerdelay=0;
triggerpulsewidth=Setup.DMD.illuminatetime+triggerdelay;
vddelay=0;
Setup.DMD.alp_returnvalue=calllib('DMD','AlpSeqTiming',Setup.DMD.deviceid,sequenceid,...
  Setup.DMD.illuminatetime, picturetime, triggerdelay, triggerpulsewidth, vddelay);
if Setup.DMD.alp_returnvalue~=0
    disp('Error setting sequence timing!');
    Setup.DMD.alp_returnvalue=0;
end

%% check current projection mode and start projection
uservar = int32(0);
uservarptr = libpointer('int32Ptr', uservar);
Setup.DMD.inquiretype=int32(2300);
[Setup.DMD.alp_returnvalue, ProjMode] = calllib('DMD', 'AlpProjInquire', ...
    Setup.DMD.deviceid, Setup.DMD.inquiretype, uservarptr);
if Setup.DMD.alp_returnvalue~=0
    disp('Error inquire sequence!');
    Setup.DMD.alp_returnvalue=0;
else
    if ProjMode==2302 %slave
        Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjStart', ...
        Setup.DMD.deviceid, sequenceid);
        if Setup.DMD.alp_returnvalue~=0
            disp('Error start constant display!');
            Setup.DMD.alp_returnvalue=0;
        end
    else
        Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjStartCont', ...
        Setup.DMD.deviceid, sequenceid);
        if Setup.DMD.alp_returnvalue~=0
            disp('Error start constant display!');
            Setup.DMD.alp_returnvalue=0;
        end
    end
end
end

