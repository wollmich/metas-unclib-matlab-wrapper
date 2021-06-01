% Example Interpolation
% Michael Wollensack METAS - 18.09.2018

% Requires MATLAB R2016 or newer

%clear all;
close all;

%% Interpolation of Sine Data
x1 = [0 1 2.5 3.6 5 7 8.1 10];
y1 = sin(x1);
% Test interpolation
res1 = TestInterpolation(x1, y1, @(x) sin(x));

%% Interpolation of Distribution
x2 = -4:4;
y2 = [0 .15 1.12 2.36 2.36 1.46 .49 .06 0];
% Test interpolation
res2 = TestInterpolation(x2, y2, @(x) NaN(size(x)));

%% Interpolation of Data
x3 = -3:3; 
y3 = [-1 -1 -1 0 1 1 1];
% Test interpolation
res3 = TestInterpolation(x3, y3, @(x) NaN(size(x)));

%% Interpolation of Oscillatory Sample Function
x4 = 0:25;
y4 = besselj(1, x4);
% Test interpolation
res4 = TestInterpolation(x4, y4, @(x) besselj(1, x));

%% Subfunction Test Interpolation
function res = TestInterpolation(x, y, f)
unc = @LinProp;
res = {};
res.x = x;
% Systematic uncertainty
res.es = unc(0, 0.1);
% Random uncertainty
res.er = unc(zeros(size(x)), 0.1.^2.*eye(numel(x)));
% Function
res.y = y + res.es + res.er;
% New X
res.xx = linspace(x(1), x(end), 20.*(numel(x) - 1) + 1);
% New Y (nominal function)
res.yy = f(res.xx);
% Linear interpolation
res.yy_1a = interpolation(x, res.y, 1, res.xx);
res.yy_1b = interpolation2(x, res.y, 1, res.xx);
res.yy_1c = interp1(x, y, res.xx); % MATLAB (only values)
% Quadratic interpolation
res.yy_2a = interpolation(x, res.y, 2, res.xx);
res.yy_2b = interpolation2(x, res.y, 2, res.xx);
% Cubic interpolation
res.yy_3a = interpolation(x, res.y, 3, res.xx);
res.yy_3b = interpolation2(x, res.y, 3, res.xx);
% Spline interpolation (not-a-knot)
res.yy_4a = spline(x, res.y, res.xx, 'not-a-knot');
res.yy_4b = spline2(x, res.y, res.xx, 'not-a-knot');
res.yy_4c = spline(x, y, res.xx); % MATLAB (only values)
% Spline interpolation (natural spline)
res.yy_5a = spline(x, res.y, res.xx, 'natural spline');
res.yy_5b = spline2(x, res.y, res.xx, 'natural spline');
% Spline interpolation (1st derivative)
res.yy_6a = spline(x, res.y, res.xx, '1st derivative', 0, '1st derivative', 0);
res.yy_6b = spline2(x, res.y, res.xx, '1st derivative', 0, '1st derivative', 0);
res.yy_6c = spline(x(:), [0; y(:); 0], res.xx); % MATLAB (only values)
% Spline interpolation (2nd derivative)
% Same result as natural spline when 2nd derivatives values are set to zero
res.yy_7a = spline(x, res.y, res.xx, '2nd derivative', 0, '2nd derivative', 0);
res.yy_7b = spline2(x, res.y, res.xx, '2nd derivative', 0, '2nd derivative', 0);
% PCHIP (not implemented in METAS UncLib, only MATLAB)
res.yy_8c = pchip(x, y, res.xx); % MATLAB (only values)
% Residuals (METAS UncLib Interpolation vs MATLAB Interpolation)
res.mres_1 = Residual(res.yy_1a, res.yy_1c);
res.mres_4 = Residual(res.yy_4a, res.yy_4c);
res.mres_6 = Residual(res.yy_6a, res.yy_6c);
% Residuals (Interpolation vs Function)
res.res_1 = Residual(res.yy_1a, res.yy);
res.res_2 = Residual(res.yy_2a, res.yy);
res.res_3 = Residual(res.yy_3a, res.yy);
res.res_4 = Residual(res.yy_4a, res.yy);
res.res_5 = Residual(res.yy_5a, res.yy);
res.res_6 = Residual(res.yy_6a, res.yy);
res.res_7 = Residual(res.yy_7a, res.yy);
res.res_8 = Residual(res.yy_8c, res.yy);
% Names
res.name_1 = 'Linear';
res.name_2 = 'Quadratic';
res.name_3 = 'Cubic';
res.name_4 = 'Spline (not-a-knot)';
res.name_5 = 'Spline (natural spline)';
res.name_6 = 'Spline (1st derivative)';
res.name_7 = 'Spline (2nd derivative)';
res.name_8 = 'PCHIP';
% Figure
res.h = figure();
subset = [1 2 4 5 6];
names = {res.name_1, res.name_2, res.name_3, res.name_4, res.name_5, res.name_6, res.name_7, res.name_8};
yy_a = [res.yy_1a; res.yy_2a; res.yy_3a; res.yy_4a; res.yy_5a; res.yy_6a; res.yy_7a; res.yy_8c];
yy_b = [res.yy_1b; res.yy_2b; res.yy_3b; res.yy_4b; res.yy_5b; res.yy_6b; res.yy_7b; res.yy_8c];

l1 = {'Points', names{subset}, 'Location', 'SouthEast'};
subplot(2,2,1);
plot(x, get_value(res.y), 'o', res.xx, get_value(yy_a(subset,:)));
title('Values of interpolation and spline');
res.a(1) = gca;
%legend(l1{:});

subplot(2,2,2);
plot(x, get_value(res.y), 'o', res.xx, get_value(yy_b(subset,:)));
title('Values of interpolation2 and spline2');
res.a(2) = gca;
%legend(l1{:});

subplot(2,2,3);
plot(x, get_stdunc(res.y), 'o', res.xx, get_stdunc(yy_a(subset,:)));
title('Uncertainties of interpolation and spline');
res.a(3) = gca;
%legend(l1{:});

subplot(2,2,4);
plot(x, get_stdunc(res.y), 'o', res.xx, get_stdunc(yy_b(subset,:)));
title('Uncertainties of interpolation2 and spline2');
res.a(4) = gca;
res.a(4).YLim = res.a(3).YLim;
res.a(3).YLim = res.a(3).YLim;
legend(l1{:});
end

function r = Residual(a ,b)
d = double(a(:)) - double(b(:));
r = d'*d;
end
