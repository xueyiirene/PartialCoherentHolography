function [] = function_feed_DMD( frames )
[LX,LY,LF] = size(frames);
disp('Loading frames on DMD')
for j = 1:(LF)
%here pick  frame for bufferload 
frame = squeeze(frames(:,:,j));
DMDvec=im2DMD(frame');
%disp(['Loading frame ' int2str(j)]);
cacheposition=j;
calllib('DMD','DLP_Img_DownloadBitplanePatternToExtMem',DMDvec(:),98304,cacheposition);
end
disp([int2str(LF) ' frames loaded'])
sequence = linspace(1,LF,LF);
calllib('DMD','DLP_Display_DisplayPatternManualForceFirstPattern');
calllib('DMD','DLP_RegIO_WriteImageOrderLut',1, sequence, numel(sequence));
calllib('DMD','DLP_Display_DisplayPatternAutoStepRepeatForMultiplePasses');
calllib('DMD','DLP_Source_SetDataSource','SL_EXT3P3');
end

