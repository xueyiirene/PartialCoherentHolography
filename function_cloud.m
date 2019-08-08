function [ frames ] = function_cloud( Setup )

LF = 64;
frames = zeros(Setup.DMD.LX,Setup.DMD.LY,LF);

for i = 1:30
x = floor(rand*(Setup.DMD.LX));
y = floor(rand*(Setup.DMD.LY));
z = 100*rand-50;
s = 20;
frames  = max(frames, function_makespot( Setup.DMD.LX,Setup.DMD.LY,x,y,z,s ));
end

f = figure(1);
for j = 1:(LF)
%here pick  frame for bufferload 
frame = squeeze(frames(:,:,j));
imagesc(frame);
pause(0.01)
end
close(f)

end

