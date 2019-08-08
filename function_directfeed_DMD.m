function [] = function_directfeed_DMD( Setup,frames )
function_feed_DMD( frames );
outputSingleScan(Setup.Daq,[0 0 0 0 0 0]);
outputSingleScan(Setup.Daq,[0 0 0 0 0 1]);
outputSingleScan(Setup.Daq,[0 0 0 0 0 0]);

end

