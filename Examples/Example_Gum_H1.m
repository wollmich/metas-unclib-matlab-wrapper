% Example Gum H1 - End-gauge calibration
% Michael Wollensack METAS - 17.12.2008

clear all;
close all;

unc = @LinProp;

disp(sprintf('\nExample GUM H1 - End-gauge calibration\n'))

%% Calibration of standard end gauge
l_s = unc(50.000623e-3, 25e-9, 1/18);

%% Measured difference between end gauges
% repeated observations
d1 = unc(215e-9, 5.8e-9, 1/24);
% random effects of comporator
d2 = unc(0, 3.9e-9, 1/5);
% systematic effects of comporator
d3 = unc(0, 6.7e-9, 1/8);
d = d1 + d2 + d3;
disp(sprintf('d      = %0.6e m', get_value(d)))
disp(sprintf('u(d)   = %0.6e m', get_stdunc(d)))
disp(sprintf('dof(d) = %0.2f', 1/get_idof(d)))

%% Thermal expansion coefficient of standard end gauge (uniform)
alpha_s = unc(11.5e-6, 2e-6/sqrt(3));

%% Temperature of test bed
% mean temperature of bed
theta_1 = unc(-0.1, 0.2);
% cyclic variation of temperature of room (arcsine)
theta_2 = unc(0, 0.5/sqrt(2));
theta = theta_1 + theta_2;

%% Difference in expansion coefficients of end gauges (uniform)
delta_alpha = unc(0, 1e-6/sqrt(3), 1/50);

%% Difference in temperatures of end gauges (uniform)
delta_theta = unc(0, 0.05/sqrt(3), 1/2);

%% Mathematical model 1
alpha = delta_alpha + alpha_s;
theta_s = theta - delta_theta;
l1 = (l_s.*(1 + alpha_s.*theta_s) + d)./(1 + alpha.*theta);
disp('Final result:')
disp(sprintf('l1      = %0.6e m', get_value(l1)))
disp(sprintf('u(l1)   = %0.6e m', get_stdunc(l1)))
disp(sprintf('dof(l1) = %0.2f', 1/get_idof(l1)))

%% Mathematical model 2
tmp1 = -l_s .* delta_alpha .* theta;
tmp2 = -l_s .* alpha_s .* delta_theta;
l2 = l_s + d + tmp1 + tmp2;
disp('Final result:')
disp(sprintf('l2      = %0.6e m', get_value(l2)))
disp(sprintf('u(l2)   = %0.6e m', get_stdunc(l2)))
disp(sprintf('dof(l2) = %0.2f', 1/get_idof(l2)))

%% Other
%get_correlation([l1 l2])
%get_unc_component(l1, [l_s d alpha_s theta delta_alpha delta_theta])
