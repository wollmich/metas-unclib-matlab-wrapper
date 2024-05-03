% Example 1 B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 07.12.2023

clear all;
close all;

unc = @LinProp;
addpath('b_least');

fprintf('\nExample B_LEAST 1\n\n')

%% Inputs
% Data files
[py, px] = b_read_cal_data('data\b_least_1_data_cal.txt');
pym = b_read_meas_data('data\b_least_1_data_meas.txt');

% Fit function
b_func = @b_linear_func;

%% Calibration
tic
[b, fweighted, res] = b_least(py, px, b_func);
toc

b_disp_cal_results(b, fweighted, res);

%% Error correction
pxm = b_eval(pym, b, b_func);

b_disp_meas_results(pxm, pym);
b_plot(py, px, pym, pxm, b, b_func);
