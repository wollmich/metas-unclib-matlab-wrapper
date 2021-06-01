% Example Gum H4 - Measurement of activity
% Michael Wollensack METAS - 17.12.2008

clear all;
close all;

unc = @LinProp;

disp(sprintf('\nExample Gum H4 - Measurement of activity\n'))

%% Counting data for determining the activity concentration
T0 = 60; % min
lambda = 1.25894e-4; % min^-1
ts = [ 243.74  984.53 1723.87 2463.17 3217.56 3956.83];
Cs = [  15380   14978   14394   13254   12516   11058];
tb = [ 305.56 1046.10 1785.43 2524.73 3279.12 4018.38];
Cb = [   4054    3922    4200    3830    3956    3980];
tx = [ 367.37 1107.66 1846.99 2586.28 3340.68 4079.94];
Cx = [  41432   38706   35860   32238   29640   26356];

Rx = (Cx - Cb)./T0.*exp(lambda.*tx);
Rs = (Cs - Cb)./T0.*exp(lambda.*ts);

data = [Rx; Rs]';

%% Analysis of data
input_values = mean(data, 1);
input_covar = cov(data)./size(data, 1);

inputs = unc(input_values, input_covar);

Rx_m = inputs(1);
Rs_m = inputs(2);
corr = get_correlation([Rx_m Rs_m]);

disp(sprintf('Rx         = %10.6f min^-1', get_value(Rx_m)))
disp(sprintf('u(Rx)      = %10.6f min^-1', get_stdunc(Rx_m)))
disp(sprintf('Rs         = %10.6f min^-1', get_value(Rs_m)))
disp(sprintf('u(Rs)      = %10.6f min^-1', get_stdunc(Rs_m)))
disp(sprintf('r(Rx,Rs)   = %10.6f', corr(1,2)))

R_m = Rx_m./Rs_m;
disp(sprintf('R          = %10.6f', get_value(R_m)))
disp(sprintf('u(R)       = %10.6f', get_stdunc(R_m)))

%% Calculation of final results
As = unc(0.1368, 0.0018); % Bq/g
ms = unc(5.0192, 0.0050); % g
mx = unc(5.0571, 0.0010); % g

Ax = As.*ms./mx*R_m;
disp('Final result:')
disp(sprintf('Ax         = %10.6f Bq/g', get_value(Ax)))
disp(sprintf('u(Ax)      = %10.6f Bq/g', get_stdunc(Ax)))
disp(sprintf('u(Ax)/Ax   = %10.6f Bq/g', get_stdunc(Ax)./get_value(Ax)))
