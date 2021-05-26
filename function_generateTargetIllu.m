function [FinalFrames] = function_generateTargetIllu(I,Mask,t1,threshold)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Mask_ref = imref2d(size(Mask)); %relate intrinsic and world coordinates
I_registered = imwarp(I,t1,'OutputView',Mask_ref);

Mask_temp=I_registered;
Mask_temp(Mask_temp<threshold)=0;
Mask_temp(Mask_temp>=threshold)=1;
FinalFrames=uint8(Mask_temp)*255;
end

