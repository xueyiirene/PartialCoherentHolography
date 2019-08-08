function [  ] = function_previewCam( cam,parameters )
h=figure;
try
while isvalid(h) 
[ Data ] = function_GetFrameCam( cam,parameters,1 );
imagesc(Data(:,:,1));  colormap gray; caxis([0 255]);
drawnow
end
end
end

