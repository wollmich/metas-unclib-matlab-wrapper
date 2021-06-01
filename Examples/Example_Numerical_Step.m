% Example Numerical Step
% Michael Wollensack METAS - 28.02.2013

clear all;
close all;

% Uncertainty inputs
a = LinProp(3, 0.3);
b = LinProp(4, 0.4);
% Compute result using linear uncertainty propagation
c = sqrt(a.*a + b.*b);
d = a + b + c;
% Compute result using a numerical function Y=f(X)
temp = numerical_step(@Example_Numerical_Function, [a b], [0.001 0.001]);
c2 = temp(1);
d2 = temp(2);
% Compare the results
jacobi = get_jacobi2([c c2 d d2], [a b])
