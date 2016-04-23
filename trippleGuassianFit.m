function [fitresult, gof] = trippleGaussianFit(x, y, z, x0, y0)
%% Customized the auto generated Curve Fitting function.
[xData, yData, zData] = prepareSurfaceData( x, y, z );

% Set up fittype and options.
ft = fittype( '(a1/sqrt(2*pi).*exp(-((x-b1).^2/2)-((y-c1).^2/2))) + (a2/sqrt(2*pi).*exp(-((x-b2).^2/2)-((y-c2).^2/2)))+(a3/sqrt(2*pi).*exp(-((x-b3).^2/2)-((y-c3).^2/2)))', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
b1 = x0-7;
b2 = x0;
b3 = x0+7;

c1 = y0+7; 
c2 = y0;
c3 = y0-7;
opts.StartPoint = [0.178132454400338 0.128014399720173 0.190433267179954 b1 b2 b3 c1 c2 c3];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, [xData, yData], zData );
% legend( h, 'untitled fit 1', 'z vs. x, y', 'Location', 'NorthEast' );
% % Label axes
% xlabel x
% ylabel y
% zlabel z
% grid on
% view( -2.3, 8.4 );


