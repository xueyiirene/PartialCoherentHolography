function [] = function_directfeed_DMD( Setup,frames )
function_feed_DMD( Setup,frames);
outputSingleScan(Setup.Daq,[0 0 0 0 0 0 1 0]);
outputSingleScan(Setup.Daq,[0 0 0 0 0 1 1 0]);
outputSingleScan(Setup.Daq,[0 0 0 0 0 0 1 0]);
end

