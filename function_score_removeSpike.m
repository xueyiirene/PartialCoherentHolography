function [ score, odata, spikeremoveflag] = function_score_removeSpike( Setup,Stim,data,spikeremoveflag)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Fs = 1/(Stim.UT(2)-Stim.UT(1));
Setup.Scorepikes.SpikeSTDDetectionthreshold=10;
Setup.RCCutoffFrequencyHz=500;
Setup.Scorepikes.MinPeakDistance=0.003;
Setup.Scorepikes.baselineduration=0.04;
Setup.Scorepikes.Photocurrentwindowduration=0.01;

% Case where we count the number of spikes... returns High pass filtered
% data and number of spikes that are away from data std.
% figure(1);plot(data);pause(0.5)
    [ HPData, ~ ] = function_RC_filter( data, 1/(Stim.UT(2)-Stim.UT(1)),Setup.RCCutoffFrequencyHz );
    data_HP = abs(HPData-mean(HPData));
    threshold = Setup.Scorepikes.SpikeSTDDetectionthreshold *std(data_HP);
    
    odata = data;
    baseline = Stim.baseline;
    data = data-baseline;
    
    nonzerovalues = find(Stim.Output(:,2));%change to (:,1) if stim light is blue, (:,2) for yellow stim
if numel( nonzerovalues)==0
    score = 0;
else

    if (~isempty(find(data_HP>threshold, 1)))
        [~,spiketimes]=findpeaks(data_HP, 'minpeakheight',threshold,'minpeakdistance',Fs*Setup.Scorepikes.MinPeakDistance);
            if ~isempty(spiketimes)
                spikeremoveflag=1;
                oridata=data;
                for k=1:length(spiketimes)
%                     [~,indspike]=min(abs(cropbegin-spiketimes(k)));
%                     cropbegin(indspike)=[];
                    if spiketimes<100
                        break;
                    else
                        if spiketimes(k)<800
                            temp=data(1:spiketimes(k)+750);
                            data(1:spiketimes(k)+750)=medfilt1(temp,500);
                        else
                            temp=data(spiketimes(k)-750:spiketimes(k)+750);
                            data(spiketimes(k)-750:spiketimes(k)+750)=medfilt1(temp,500);
                        end
                        end
                end
            end
    end
  
        score = max(medfilt1(data));
end
end

