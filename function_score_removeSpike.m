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

    [ HPData, ~ ] = function_RC_filter( data, 1/(Stim.UT(2)-Stim.UT(1)),Setup.RCCutoffFrequencyHz );
%     figure(1);plot(data);pause(0.2);
    data_HP = abs(HPData-mean(HPData));
    threshold = Setup.Scorepikes.SpikeSTDDetectionthreshold *std(data_HP);
    
    odata = data;
    baseline = mean(data(1:floor(Fs*Setup.Scorepikes.baselineduration)));
    data = data-baseline;
    baselinepeak = max(abs(medfilt1(data(1:floor(Fs*Setup.Scorepikes.baselineduration)))));
   
    nonzerovalues = find(Stim.Output(:,1));
if numel( nonzerovalues)==0
    score = 0;
else
    cropbeginInd=1:length(nonzerovalues)/Stim.Npeaks:length(nonzerovalues);
    cropbegin = nonzerovalues(cropbeginInd);
    if (~isempty(find(data_HP(nonzerovalues(1):nonzerovalues(end))>threshold, 1)))
        [~,spiketimes]=findpeaks(data_HP, 'minpeakheight',threshold,'minpeakdistance',Fs*Setup.Scorepikes.MinPeakDistance);
            if ~isempty(spiketimes)
                spikeremoveflag=1;
                figure();plot(data);
                for k=1:length(spiketimes)
                    [~,indspike]=min(abs(cropbegin-spiketimes(k)));
                    cropbegin(indspike)=[];
                end
            end
    end
    cropend = cropbegin+ floor(Fs*Setup.Scorepikes.Photocurrentwindowduration);
    score_temp=zeros(size(cropbegin));
    for k=1:length(cropbegin)
        cropdata = data(cropbegin(k):cropend(k));
        % blue light = Stim.Output(:,1)
        %%%%%%%Setup.Scorepikes.Photocurrentwindowduration
    %     figure(1);subplot(1,2,1);cla;plot(data);hold on;plot(Stim.Output(:,1)/100);
    %     subplot(1,2,2);plot(cropdata);
    %     pause(1)
        score_temp(k) = max(abs(medfilt1(cropdata))-baselinepeak);
    end
    score=mean(score_temp);
end
end

