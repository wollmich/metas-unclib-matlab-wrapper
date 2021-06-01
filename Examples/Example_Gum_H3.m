% Example Gum H3 - Calibration of a thermometer
% Michael Wollensack METAS - 17.12.2008

clear all;
close all;

unc = @LinProp;

disp(sprintf('\nExample Gum H3 - Calibration of a thermometer\n'))

%% Data used to obtain a linear calibration curve for a thermometer
t0 = 20;
% Thermometer reading
tk = [21.521 22.012 22.512 23.003 23.507 23.999 24.513 25.002 25.503 26.010 26.511];
% Observed correction
bk = [-0.171 -0.169 -0.166 -0.159 -0.164 -0.165 -0.156 -0.157 -0.159 -0.161 -0.160];


%% Least-square fitting (see MATLAB help polyfit)
[p s] = polyfit((tk - t0), bk, 1);
Rinv = inv(s.R);
cv = (Rinv*Rinv')*s.normr^2/s.df;

tmp = unc(p, cv);
y1 = tmp(2);
y2 = tmp(1);
corr = get_correlation([y1 y2]);

disp(sprintf('y1         = %10.6f °C', get_value(y1)))
disp(sprintf('u(y1)      = %10.6f °C', get_stdunc(y1)))
disp(sprintf('y2         = %10.6f', get_value(y2)))
disp(sprintf('u(y2)      = %10.6f', get_stdunc(y2)))
disp(sprintf('r(y1,y2)   = %10.6f', corr(1,2)))

%% Predicted correction and difference between observed and
b_tk = y1 + y2.*(tk - t0);
delta_b = bk - b_tk;

%% Uncertainty of predicted value (Example t = 30°C)
b_30 = y1 + y2.*(30 - t0);
disp(sprintf('b(30°C)    = %10.6f °C', get_value(b_30)))
disp(sprintf('u(b(30°C)) = %10.6f °C', get_stdunc(b_30)))

%% Uncertainty of predicted value (Plot t = 20...30°C)
t = 20:0.01:30;
b = y1 + y2.*(t - t0);

figure();
title('Example Gum H3 - Calibration of a thermometer');
xlabel('t (°C)');
ylabel('b (°C)');
grid on;
hold on;
value = get_value(b);
unc = get_stdunc(b);
v1 = value + unc;
v2 = value - unc;
fill([t,t(end:-1:1)],[v1,v2(end:-1:1)]','m');
plot(t,value,'b-');

plot(tk, bk, 'bo');
hold off
