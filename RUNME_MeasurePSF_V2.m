clear all;close all;clc
%%
[ Setup ] = Function_Load_Parameters();
%%
Result.Zrange = 200; % In microns, range of Zstack
Result.Nsteps = 20; % Number of slices in Z stack
ProjMode = function_CheckProjMode_DMD(Setup);
if ProjMode==2302
    Setup=function_StopProj_DMD(Setup);
    Setup= function_DMDProjMode(Setup,'master');
end
brightfield = uint8(ones(1600, 2560)*255);
function_directfeed_DMD( Setup,brightfield);
outputSingleScan(Setup.Daq,[1.2 0 0 0 0 0 1]);
%% control stage
xyz_um = getPosition(Setup.SutterStage);
ZDepths = linspace(0,Result.Zrange,Result.Nsteps);
for zz = 1:Result.Nsteps
    moveTo(Setup.SutterStage,[xyz_um(1);xyz_um(2);xyz_um(3)+ZDepths(zz)]);
    pause(1);
    outputSingleScan(Setup.Daq,[2,0,0,0,0,0,0]);
    outputSingleScan(Setup.Daq,[2,0,0,0,0,1,0]);
%     outputSingleScan(Setup.Daq,[3/200*zz+1,0,0,0,0,1,0]);
    pause(1);
end
outputSingleScan(Setup.Daq,[0,0,0,0,0,0,0]);
%% back to original location
moveTo(Setup.SutterStage,[xyz_um(1);xyz_um(2);xyz_um(3)]);
%% Stop projection
[Setup]=function_StopProj_DMD(Setup);
function_Stop_DMD(Setup);
%% Stop piezo
fclose(Setup.SutterStage);