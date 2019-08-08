function [ score, odata] = function_scorespikes( Setup,Stim,data)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Fs = 1/(Stim.UT(2)-Stim.UT(1));

% Case where we count the number of spikes... returns High pass filtered
% data and number of spikes that are away from data std.
if Setup.Scorepikes.Method == 1
    [ HPData, ~ ] = function_RC_filter( data, 1/(Stim.UT(2)-Stim.UT(1)),Setup.RCCutoffFrequencyHz );
    data = HPData-mean(HPData);
    odata = data;
    data = abs(data);
    threshold = Setup.Scorepikes.SpikeSTDDetectionthreshold *std(data);
    if (~isempty(find(data>threshold, 1)))
        [~,spiketimes]=findpeaks(data, 'minpeakheight',threshold,'minpeakdistance',Fs*Setup.Scorepikes.MinPeakDistance);
        spiketimes = spiketimes';
        spiketimes = spiketimes/Fs;
    else
        spiketimes = [];
    end
    spikes1=find(spiketimes>0);
    score=length(spikes1);
    
    %Case where we measure the magnitude of the biggest peak.
elseif Setup.Scorepikes.Method == 0
    odata = data;
    
    baseline = mean(data(1:floor(Fs*Setup.Scorepikes.baselineduration)));
    data = data-baseline;
    nonzerovalues = find(Stim.Output(:,1));
if numel( nonzerovalues)==0
    score = 0;
else
    cropbegin = nonzerovalues(1);
    cropend = cropbegin+ floor(Fs*Setup.Scorepikes.Photocurrentwindowduration);
    cropdata = data(cropbegin:cropend);
    % blue light = Stim.Output(:,1)
    %%%%%%%Setup.Scorepikes.Photocurrentwindowduration
    score = max(abs(cropdata)); 
end

else
    %Case where we are being stupid
    disp('Not a valid spike evaluation method')
end
end

