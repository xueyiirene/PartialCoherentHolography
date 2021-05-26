function [Setup]=function_StopDMDSequence(Setup)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Setup.DMD.alp_returnvalue=0;
Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpProjHalt', Setup.DMD.deviceid);
    if Setup.DMD.alp_returnvalue~=0
        disp('Error stop projection!');
        Setup.DMD.alp_returnvalue=0;
    else
        disp('Stop projecting!');
    end
    
for n=1:numel(Setup.DMD.sequenceid)
    
    Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpSeqFree', Setup.DMD.deviceid, Setup.DMD.sequenceid(n));
    if Setup.DMD.alp_returnvalue~=0
        disp('Error stop current DMD sequences!');
        Setup.DMD.alp_returnvalue=0;
    else
        disp(['DMD sequences' num2str(Setup.DMD.sequenceid(n)) 'are removed!']);
        Setup.DMD.sequenceid(n) = [];
    end
end
end

