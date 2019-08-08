function [ d,b ] = function_find_Nearest_blob( frame,radius, xxx,yyy,D )
%find the block that is nearest from the xx,yy location (find a peak near
%estimated value, within distance D
smoothframe = imgaussfilt(frame,radius);
[LX,LY] = size(smoothframe);
UX = 1:LX;
UY = 1:LY;
[XX,YY] = ndgrid(UX,UY);
killme = double((XX-xxx).^2 + (YY-yyy).^2 >D^2);
smoothframe(killme>0) = 0;
[~,b] = max(max(smoothframe));
[~,d] = max(smoothframe(:,b));
end

