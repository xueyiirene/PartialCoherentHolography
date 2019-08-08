
function [ HPData, LPData ] = function_RC_filter( data, samplingfrequency,CutoffFrequency )
rollingaveragesize = samplingfrequency/CutoffFrequency;
LPData = smoothdata(data,'gaussian',rollingaveragesize);
HPData = data-LPData;

end