function [ Data ] = function_GetFrameCam( cam,Nico,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Data = zeros(Nico.H,Nico.W,Nico.Bits/8,n);
for i = 1:n
cam.Acquisition.Freeze(uc480.Defines.DeviceParameter.Wait);
[~,tmp] = cam.Memory.CopyToArray(Nico.MemId);
tmps = reshape(uint8(tmp),[Nico.Bits/8,Nico.W,Nico.H]);
Data(:,:,:,i) = permute(tmps,[3,2,1]);
end
Data = squeeze(Data);
end

