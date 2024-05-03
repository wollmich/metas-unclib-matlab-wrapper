% Example 1 B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 07.12.2023

clear all;
close all;

unc = @LinProp;
addpath('b_least');

fprintf('\nExample B_LEAST 1\n\n')

%% Inputs
% Calibration inputs
x1 = unc(NormalDistribution(4.5, 0.045), 'Reference gas No. 1 x');
x2 = unc(NormalDistribution(18.75, 0.1875), 'Reference gas No. 2 x');
x3 = unc(NormalDistribution(50.0, 0.5), 'Reference gas No. 3 x');

y1 = unc(NormalDistribution(0.1969, 0.003938), 'Reference gas No. 1 y');
y2 = unc(NormalDistribution(0.7874, 0.015748), 'Reference gas No. 2 y');
y3 = unc(NormalDistribution(2.0228, 0.040456), 'Reference gas No. 3 y');

px = [x1; x2; x3];
py = [y1; y2; y3];

% Measurement inputs
ym1 = unc(NormalDistribution(2.58e-1, 5.16e-3), 'Mixture No. 1 y');
ym2 = unc(NormalDistribution(6.0e-1, 1.2e-2), 'Mixture No. 2 y');
ym3 = unc(NormalDistribution(1.8e0, 3.6e-2), 'Mixture No. 3 y');

pym = [ym1; ym2; ym3];

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
