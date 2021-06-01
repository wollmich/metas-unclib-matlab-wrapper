% Example Right Triangle
% Michael Wollensack METAS - 18.12.2008

clear all;
close all;

unc = @LinProp;

disp(sprintf('\nExample Right Triangle\n'))

%% Definition of the inputs
a = unc(3.0, 0.3);
b = unc(4.0, 0.4);

%% Compute the outputs
c = sqrt(a.*a + b.*b);
U = a + b + c;
A = a.*b./2;

%% Display some results
disp(sprintf('c      = %7.3f', get_value(c)))
disp(sprintf('u(c)   = %7.3f', get_stdunc(c)))

disp(sprintf('U      = %7.3f', get_value(U)))
disp(sprintf('u(U)   = %7.3f', get_stdunc(U)))

disp(sprintf('A      = %7.3f', get_value(A)))
disp(sprintf('u(A)   = %7.3f', get_stdunc(A)))

temp = get_correlation([A U]);
disp(sprintf('r(A,U) = %7.3f', temp(1,2)))
