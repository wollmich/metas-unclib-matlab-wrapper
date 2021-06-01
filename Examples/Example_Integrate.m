% Example Integrate
% Michael Wollensack METAS - 17.09.2018

%clear all;
%close all;

unc = @LinProp;

%% Integral of sin(x) from 0 to pi/2
n1 = 5;
n2 = 101;
x1 = linspace(0, pi./2, n1);
x2 = linspace(0, pi./2, n2);

y1 = sin(x1);
y2 = sin(x2);

y1b = unc(y1, 1e-3.*eye(length(y1)));
y2b = unc(y2, 1e-3.*eye(length(y2)));

y3b = interpolation(x1, y1b, 1, x2);
y4b = interpolation(x1, y1b, 2, x2);
y5b = interpolation(x1, y1b, 3, x2);
y6b = spline(x1, y1b, x2);

% Integral of function
r = -cos(pi./2) + cos(0)
% Integral of values using cumtrapz of n1 data points
r1_0 = cumtrapz(x1, y1);
r1_0_end = r1_0(end)
% Integral (linear) of n1 data points
r1_1 = integrate(x1, y1b, 1);
r1_1_end = r1_1(end)
d1_1_1_0 = get_value(r1_1) - r1_0;
d1_1_1_0_end = d1_1_1_0(end)
% Integral (quadratic) of n1 data points
r1_2 = integrate(x1, y1b, 2);
r1_2_end = r1_2(end)
% Integral (cubic) of n1 data points
r1_3 = integrate(x1, y1b, 3);
r1_3_end = r1_3(end)
% Integral (spline) of n1 data points
r1_s = splineintegrate(x1, y1b);
r1_s_end = r1_s(end)
% Integral of values using cumtrapz of n2 data points
r2_0 = cumtrapz(x2, y2);
r2_0_end = r2_0(end)
% Integral (linear) of n2 data points
r2_1 = integrate(x2, y2b, 1);
r2_1_end = r2_1(end)
d2_1_2_0 = get_value(r2_1) - r2_0;
d2_1_2_0_end = d2_1_2_0(end)
% Integral (quadratic) of n2 data points
r2_2 = integrate(x2, y2b, 2);
r2_2_end = r2_2(end)
% Integral (cubic) of n2 data points
r2_3 = integrate(x2, y2b, 3);
r2_3_end = r2_3(end)
% Integral (spline) of n2 data points
r2_s = splineintegrate(x2, y2b);
r2_s_end = r2_s(end)
% Integral (linear) of interpolated (linear) data points
r3_1 = integrate(x2, y3b, 1);
r3_1_end = r3_1(end)
d3_1_1_1_end = r3_1(end) - r1_1(end)
% Integral (quadratic) of interpolated (quadratic) data points
r4_2 = integrate(x2, y4b, 2);
r4_2_end = r4_2(end)
d4_2_1_2_end = r4_2(end) - r1_2(end)
% Integral (cubic) of interpolated (cubic) data points
r5_3 = integrate(x2, y5b, 3);
r5_3_end = r5_3(end)
d5_3_1_3_end = r5_3(end) - r1_3(end)
% Integral (spline) of interpolated (spline) data points
r6_s = splineintegrate(x2, y6b);
r6_s_end = r6_s(end)
d6_s_1_s_end = r6_s(end) - r1_s(end)
