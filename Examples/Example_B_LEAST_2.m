% Example 2 B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 07.12.2023

clear all;
close all;

unc = @LinProp;
addpath('b_least');

fprintf('\nExample B_LEAST 2\n\n')

%% Inputs
% Data files
[py, px] = b_read_cal_data('data\b_least_2_data_cal.txt');
pym = b_read_meas_data('data\b_least_2_data_meas.txt');

% Fit function
b_func1 = @b_linear_func;
b_func2 = @b_second_order_ploy;

%% Calibration 1
fprintf('\nLinear function\n\n')
tic
[b1, fweighted1, res1] = b_least(py, px, b_func1);
toc

b_disp_cal_results(b1, fweighted1, res1);

%% Error correction 1
pxm1 = b_eval(pym, b1, b_func1);

b_disp_meas_results(pxm1, pym);
b_plot(py, px, pym, pxm1, b1, b_func1);

%% Calibration 2
fprintf('\nSecond order polynomial\n\n')
tic
[b2, fweighted2, res2] = b_least(py, px, b_func2);
toc

b_disp_cal_results(b2, fweighted2, res2);

%% Error correction 2
pxm2 = b_eval(pym, b2, b_func2);

b_disp_meas_results(pxm2, pym);
b_plot(py, px, pym, pxm2, b2, b_func2);
