% Example Roots
% Michael Wollensack METAS - 15.08.2018

%clear all;
%close all;

%% Roots of Quadratic Polynom
% 3.*x.^2 - 2.*x - 4 = 0

p2 = LinProp(3, 0.3, 'p2');
p1 = LinProp(-2, 0.2, 'p1');
p0 = LinProp(-4, 0.4, 'p0');
p = [p2 p1 p0];

r = roots(p)
r_value = roots(get_value(p))

%% Root of Quartic Polynom
% x.^4 - 1 = 0
q4 = LinProp(1, 0.1, 'q4');
q3 = LinProp(0, 0.1, 'q3');
q2 = LinProp(0, 0.1, 'q2');
q1 = LinProp(0, 0.1, 'q1');
q0 = LinProp(-1, 0.1, 'q0');
q = [q4 q3 q2 q1 q0];

s = roots(q)
s_value = roots(get_value(q))
