function [ d,b ] = function_find_blob( frame,radius )
smoothframe = imgaussfilt(frame,radius);
[~,b] = max(max(smoothframe));
[~,d] = max(smoothframe(:,b));


end

