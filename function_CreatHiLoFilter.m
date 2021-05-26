function [fw,fhp,flp]=function_CreatHiLoFilter(I,Grid_period)
%Create HiLo filter to remove scattered photons
%   Yi Xue, 2021
R=1;%originally 10. Gaussian filter constant
FilterPower=1;%super Gaussian filter =2; Gaussian filter =1
numRows=size(I,1);numCols=size(I,2);
lp_sigma    = floor(Grid_period*2);   % pixels
w_sigma     = floor(Grid_period*4);   % pixels
border      = 5*lp_sigma;             % pixels
% Create filter for HiLo,1500x1500pixel image with 19pixel grid
lp_kernel = zeros(numRows+2*border,numCols+2*border);
w_kernel  = zeros(numRows+2*border,numCols+2*border);
% In this algorithm, hp filter and w filter are duplicated but has
% different standard deviation (lp_sigma and w_sigma)
one_over_two_sigma2_lp = 1/(2*lp_sigma^2);
one_over_two_sigma2_w = 1/(2*w_sigma^2);
for rr=0:(numRows+2*border-1)
    if rr<floor((numRows+2*border)/2)
        yVal2 = rr*rr;
    else
        yVal2 = (numRows+2*border-rr)*(numRows+2*border-rr);
    end
    for cc=0:(numCols+2*border-1)
        if cc<floor((numCols+2*border)/2)
            xVal2 = cc*cc;
        else
            xVal2 = (numCols+2*border-cc)*(numCols+2*border-cc);
        end
        % Gaussian distribution filters
        lp_kernel(1+rr,1+cc)=exp(-R*((xVal2+yVal2)*one_over_two_sigma2_lp).^FilterPower);
        w_kernel(1+rr,1+cc)=exp(-R*((xVal2+yVal2)*one_over_two_sigma2_w).^FilterPower);
    end
end

lp_kernel=lp_kernel/sum(lp_kernel(:));
lp_kernel_k=fft2(lp_kernel);
hp_kernel_k=1-lp_kernel_k;
hp_kernel=ifft2(hp_kernel_k);
% lp_kernel_k_line=abs(fftshift(lp_kernel_k(1,:)));
% hp_kernel_k_line=abs(fftshift(hp_kernel_k(1,:)));

w_kernel=w_kernel/sum(w_kernel(:));
% w_kernel_k=fft2(w_kernel);
% w_kernel_k=1-w_kernel_k;
% w_kernel=ifft2(w_kernel_k);
% w_kernel_k_line=abs(fftshift(w_kernel_k(1,:)));

fw=fft2(w_kernel);
fhp=fft2(hp_kernel);
flp=fft2(lp_kernel);
% 
% fw=gpuArray(fw);
% fhp=gpuArray(fhp);
% flp=gpuArray(flp);
end

