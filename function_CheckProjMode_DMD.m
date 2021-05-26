function ProjMode = function_CheckProjMode_DMD(Setup)
uservar = int32(0);
uservarptr = libpointer('int32Ptr', uservar);
Setup.DMD.inquiretype=int32(2300);
[Setup.DMD.alp_returnvalue, ProjMode] = calllib('DMD', 'AlpProjInquire', ...
    Setup.DMD.deviceid, Setup.DMD.inquiretype, uservarptr);
if Setup.DMD.alp_returnvalue~=0
    disp('Error inquire sequence!');
    Setup.DMD.alp_returnvalue=0;
end
end

