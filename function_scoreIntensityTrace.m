function [score] = function_scoreIntensityTrace(Stim,data)
%Find the maximum value of intensity after each photo-stimulation
%   Yi Xue, 2021
% odata = data;
% k=(mean(data(end-49:end))-mean(data(1:50)))/(numel(data)-1);
% baseline=k*((1:numel(data))-1)+mean(data(1:50));
%     baseline = mean(data(1:floor(Fs*Setup.Scorepikes.baselineduration)));
% if size(data,1)~=size(baseline,1)
%     baseline=baseline';
% end
% data = data-baseline;

% if numel(Stim.nonzerovalues)==0
%     score = 0;
% else
    %%%%%%Setup.Scorepikes.Photocurrentwindowduration
    score = max(abs(data)); 
% end
end

