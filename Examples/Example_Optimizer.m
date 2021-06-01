% Example Optimizer
% Michael Wollensack METAS - 29.09.2017

clear all;
close all;

% Start values for the variable optimization parameters
xStart = [2.5, 2];
% Constant optimization parameters
p = [LinProp(6, 2), LinProp(-4, 1), LinProp(-16, 1), LinProp(2, 1)];

% Analytic solution for example objective function (just for comparison)
x0 = [(-p(2) + sqrt(p(2) * p(2) - 4 * p(1) * p(3))) / (2 * p(1));...
      log(p(4))];      

% Optimize using Metas.UncLib.Optimization (Levenberg-Marquardt)
tic;x1 = optimizer(@Example_Objective_Function, xStart, p);toc

% Optimize using Metas.UncLib.Optimization (Trust-Region)
tic;x2 = optimizer(@Example_Objective_Function, xStart, p, true, [], [], [], 0, 'TrustRegion');toc

% Optimize using MATLAB Optimization Toolbox
tic;x3 = matlab_optimizer(@Example_Objective_Function, xStart, p);toc

