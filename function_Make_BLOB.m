function [ frame ] = function_Make_BLOB( LX,LY,x,y,s )
frame = zeros(LX,LY);
dx = floor(x-s/2);
dy = floor(y-s/2);
ex = floor(x+s/2);
ey = floor(y+s/2);
dx = max(dx,1);
dy = max(dy,1);
ex = min(ex,LX);
ey = min(ey,LY);
frame(dx:ex,dy:ey) = 1;

end

