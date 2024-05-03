% Example 3 B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 07.12.2023

clear all;
close all;

unc = @LinProp;
addpath('b_least');

fprintf('\nExample B_LEAST 3\n\n')

%% Inputs
% Data files
[py, px] = b_read_cal_data('data\b_least_3_data_cal.txt');
pym = b_read_meas_data('data\b_least_3_data_meas.txt');

% Fit function
b_func1 = @b_linear_func;
b_func2 = @b_power_func;
b_func3 = @b_exp_func;

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
fprintf('\nPower function\n\n')
tic
[b2, fweighted2, res2] = b_least(py, px, b_func2);
toc

b_disp_cal_results(b2, fweighted2, res2);

%% Error correction 2
pxm2 = b_eval(pym, b2, b_func2);

b_disp_meas_results(pxm2, pym);
b_plot(py, px, pym, pxm2, b2, b_func2);

%% Calibration 3
fprintf('\nExponential function\n\n')
tic
[b3, fweighted3, res3] = b_least(py, px, b_func3);
toc

b_disp_cal_results(b3, fweighted3, res3);

%% Error correction 3
pxm3 = b_eval(pym, b3, b_func3);

b_disp_meas_results(pxm3, pym);
b_plot(py, px, pym, pxm3, b3, b_func3);
