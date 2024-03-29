function [fitresult, gof] = createFit(UX, temp)
%CREATEFIT(UX,TEMP)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : UX
%      Y Output: temp
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 10-Jul-2019 17:27:18


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( UX, temp );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0];
opts.StartPoint = [0.250559042216675 0 25.61677657228];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'temp vs. UX', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel UX
% ylabel temp
% grid on


