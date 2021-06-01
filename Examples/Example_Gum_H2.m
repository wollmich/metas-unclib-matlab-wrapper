% Example Gum H2 - Simultaneous resistance and reactance measurement
% Michael Wollensack METAS - 18.12.2008

clear all;
close all;

unc = @LinProp;

disp(sprintf('\nExample Gum H2 - Simultaneous resistance and reactance measurement\n'))

%% Definition of the inputs

meas = [5.007 19.663e-3 1.0456; ...
        4.994 19.639e-3 1.0438; ...
        5.005 19.640e-3 1.0468; ...
        4.990 19.685e-3 1.0428; ...
        4.999 19.678e-3 1.0433];

input_values = mean(meas, 1);
input_covar = cov(meas)./size(meas, 1);

inputs = unc(input_values, input_covar);

v = inputs(1);
i = inputs(2);
phi = inputs(3);

%% Compute the outputs

r = v./i.*cos(phi);
x = v./i.*sin(phi);
z = v./i;

outputs = [r x z];

output_values = get_value(outputs)
output_stdunc = get_stdunc(outputs)
output_corr = get_correlation(outputs)
