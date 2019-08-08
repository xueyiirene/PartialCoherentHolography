function [ cam,Nico ] = function_StartCam( )
NET.addAssembly('C:\Program Files\Thorlabs\Scientific Imaging\DCx Camera Support\Develop\DotNet\uc480DotNet.dll');
cam = uc480.Camera;
cam.Init(0);
cam.Display.Mode.Set(uc480.Defines.DisplayMode.DiB);
cam.Trigger.Set(uc480.Defines.TriggerMode.Software);
[~, Nico.MemId] = cam.Memory.Allocate(true);
[~, Nico.W,Nico.H,Nico.Bits,~] = cam.Memory.Inquire(Nico.MemId);
end

