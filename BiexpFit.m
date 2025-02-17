function [fitresult, gof, output] = BiexpFit(b, data, op)
%CREATEFIT(B,MEDULLA_1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : b
%      Y Output: Medulla_1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
% 
%  fitting Parameters 
%  d = S_0
%  a = f_fast
%  b = D_slow
%  c = D_fast

%  Auto-generated by MATLAB on 19-Aug-2020 15:23:59


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( b, data );

% Set up fittype and options.
ft = fittype( 'd*((1-a)*exp(-b*x) + a*exp(-c*x))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%opts.Algorithm = 'Levenberg-Marquardt';
%opts.Algorithm = 'Trust-Region';
%opts.DiffMinChange = 1e-07;
opts.Display = 'Off';

opts.Algorithm = op.Algorithm;
opts.Display = op.Display;
opts.StartPoint = op.StartPoint ;

if strcmp(opts.Algorithm, 'Trust-Region')
    opts.Lower =      op.Lower;
    opts.Upper =      op.Upper;
end

% opts.Lower = [0.001 0.001 0.00001 0.002 0.009];
% opts.MaxIter = 500;
% opts.StartPoint = [0.16 0.06 0.0019 0.0097 0.551];
% opts.Upper = [0.999 0.999 0.007 0.05 1];

% opts.Lower = [ 0.001 0.001 0.0007 0.0025 0.015];
% opts.MaxIter = 500;
% opts.StartPoint = [ 0.4 0.03 0.0009 0.0058 0.16];
% opts.Upper = [ 0.999 0.999 0.0025 0.015 0.2];

% Fit model to data.
[fitresult, gof, output] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'data vs. b', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel b
% ylabel data
% grid on


