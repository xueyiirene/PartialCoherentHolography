function [ FWHM ] = function_FWHMofGaussian( Gfit, X )
%Calcualte FWHM from Gaussian fit
%   Gfit: cfit object
x=linspace(X(1),X(end),length(X)*10);
y=Gfit.a1*exp(-((x-Gfit.b1)/Gfit.c1).^2);
[peak,ind]=max(y);
[~,ind1]=min(abs(y(1:ind)-peak/2));
[~,ind2]=min(abs(y(ind:end)-peak/2));
FWHM=x(ind2+ind-1)-x(ind1);
end

