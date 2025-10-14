% Example ODE Solver
% Michael Wollensack METAS - 14.10.2025

clear all;
close all;

% Uncertainty inputs
y = [LinProp(1)];
x = [LinProp(0, 0.1); LinProp(1, 0.1); LinProp(2, 0.1); LinProp(3, 0.1)];
p = [LinProp(-1.2, 0.01)];
% Compute result 
yres = exp(p(1)*x)
% Compute result using ODE solver
ytbl = ode_solver(y, x, p, 1e-12, @Example_ODE_Function)
