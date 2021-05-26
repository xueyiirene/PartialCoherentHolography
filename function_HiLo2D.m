function [HiLo] = function_HiLo2D(Uniform_image, Structured_image, Grid_period,fw,fhp,flp)
Iu_pre = Uniform_image;
Is_pre = Structured_image;
[numRows, numCols] = size(Iu_pre);

grid_period = Grid_period; 
lp_sigma    = floor(grid_period*2);   % pixels
% w_sigma     = floor(grid_period*4);   % pixels
border      = 5*lp_sigma;             % pixels
eta         = 0.5;                    % original : 4   
remove_noise = 1;                     % boolean
%% step 2: take care of borders

if(remove_noise)    % remove salt and pepper noise
    Iu_pre = medfilt2(Iu_pre,[3 1]);
    Is_pre = medfilt2(Is_pre,[3 1]);
end

Iu_raw=zeros(numRows+2*border,numCols+2*border);
Is_raw=zeros(numRows+2*border,numCols+2*border);
Iu_raw(1+border:end-border,1+border:end-border)=double(Iu_pre);
Is_raw(1+border:end-border,1+border:end-border)=double(Is_pre);
%% step 3: hilo algorithm
meanIu = mean2(Iu_raw(1+border:end-border,1+border:end-border));
meanIs = mean2(Is_raw(1+border:end-border,1+border:end-border));
Iu_n  = Iu_raw/meanIu;
Is_n  = Is_raw/meanIs;
Id    = (Is_n-Iu_n)*pi/2;
% Idw   = ifft2(fft2(Id).*fw);
% Idabs = abs(Idw);

% run on gpu
Id=gpuArray(Id);
% meanIu=gpuArray(meanIu);
Iu_raw=gpuArray(Iu_raw);

% Ilow  = ifft2(fft2(Id).*fw.*flp);
Ilow  = ifft2(fft2(Id).*fw.*flp)*meanIu;%bandpass filter can remove very low 
%frequency background fluorescence, similar to background substruction. 
% Ilow  = ifft2(fft2(Id).*flp)*meanIu;
Ihigh = ifft2(fft2(Iu_raw).*fhp);
% %normalize Ilow and Ihigh
% Ilow=Ilow/max(max(abs(Ilow)));
% Ihigh=Ihigh/max(max(abs(Ihigh)));

Ihigh=gather(Ihigh);%copy Ihigh to matlab workspace
Ilow=gather(Ilow);

% Ilow(Ilow>0)=0;

Ihilo = Ihigh + eta*Ilow;
%remove real part of Ihilo smaller than zero. This line removes the
%background around objects.
Ihiloz = Ihilo;
Ihiloz(Ihilo<0) = 0;
HiLo = Ihiloz(1+border:border+numRows,1+border:border+numCols);
