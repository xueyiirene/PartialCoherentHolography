function [Setup] = function_StartProj_DMD(Setup, sequenceid)
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

