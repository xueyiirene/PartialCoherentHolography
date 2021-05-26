function [Setup] = function_DMDProjMode(Setup,ProjTrigMode)
switch ProjTrigMode
    case 'master'
        Setup.DMD.TriggerModeValue=2301;
    case 'slave'
        Setup.DMD.TriggerModeValue=2302;
end
    Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjControl', ...
    Setup.DMD.deviceid, Setup.DMD.TriggerMode, Setup.DMD.TriggerModeValue);
    if Setup.DMD.alp_returnvalue~=0
        if Setup.DMD.TriggerModeValue==2302
            disp('Error running DMD as slave!');
        elseif Setup.DMD.TriggerModeValue==2301
            disp('Error running DMD as master!');
        end
        Setup.DMD.alp_returnvalue=0;
    else
        uservar = int32(0);
        uservarptr = libpointer('int32Ptr', uservar);
        Setup.DMD.inquiretype=int32(2300);
        [Setup.DMD.alp_returnvalue, ProjMode] = calllib('DMD', 'AlpProjInquire', ...
            Setup.DMD.deviceid, Setup.DMD.inquiretype, uservarptr);
        if Setup.DMD.alp_returnvalue~=0
            disp('Error inquire sequence!');
            Setup.DMD.alp_returnvalue=0;
        else
            if ProjMode==2302
                disp('Running DMD as a slave!');
            elseif ProjMode==2301
                disp('Running DMD as a master!');
            end
        end
    end
end

