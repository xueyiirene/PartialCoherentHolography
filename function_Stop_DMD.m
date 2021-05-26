function [] = function_Stop_DMD(Setup)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Setup.DMD.alp_returnvalue=0;
Setup.DMD.alp_returnvalue = calllib('DMD', 'AlpDevFree', Setup.DMD.deviceid);
if Setup.DMD.alp_returnvalue~=0
    disp('Error to stop DMD!');
else
    disp('DMD is stopped!');
end

end

