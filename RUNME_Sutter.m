clear all;close all;clc
%%
obj = sutterMP285('COM5');
xyz_um = getPosition(obj);
%%
dx=1;dy=1;dz=1;
moveTime = moveTo(obj,[xyz_um(1)+dx;xyz_um(2)+dy;xyz_um(3)+dz]);
%%
fclose(obj);