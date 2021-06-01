% Example Eigenvalue Problems
% Michael Wollensack METAS - 14.08.2018

%clear all;
%close all;

%% Symmetric Eigenvalue problem
s = LinProp(0, 0.1);
e11 = LinProp(1, 0.1);
e21 = LinProp(2, 0.1);
e22 = LinProp(5, 0.1);

m1 = [e11 e21;e21 e22] + s;

[v1 d1] = eig(m1)
[v1_values d1_value] = eig(get_value(m1))
% Check
check1 = m1*v1 - v1*d1

%% Non-symmetric Eigenvalue Problem
m2 = [e11 0.1.*e21;e21 e22] + s;

[v2 d2] = eig(m2)
[v2_values d2_value] = eig(get_value(m2))
% Check
check2 = m2*v2 - v2*d2

%% Linear Eigenvalue Problem
a0 = [e11 e21;e21 e22] + s;
a1 = [1 2;0 3] + s;

[v3 d3] = eig(a0, a1)
% Check
check3 = a0*v3 + a1*v3*d3

%% Quadratic Eigenvalue Problem
a2 = [4 0;5 6] + s;

[v4 d4] = eig(a0, a1, a2)
% Check
check4 = a0*v4 + a1*v4*d4 + a2*v4*d4.^2

%% Cubic Eigenvalue Problem
a3 = [7 8;9 0] + s;

[v5 d5] = eig(a0, a1, a2, a3)
% Check
check5 = a0*v5 + a1*v5*d5 + a2*v5*d5.^2 + a3*v5*d5.^3

%% Over-determined Linear Eigenvalue Problem
b0 = [e11 e21;e21 e22;2 3] + s;
b1 = [4 5;6 7;8 9] + s;

[v6 d6] = eig(b0, b1)
% Check
check6 = b0'*b0*v6 + (b0'*b1 + b1'*b0)*v6*d6 + b1'*b1*v6*d6.^2

%% Over-determined Quadratic Eigenvalue Problem
b2 = [10 11;12 13;14 15] + s;

[v7 d7] = eig(b0, b1, b2)
% Check
check7 = b0'*b0*v7 + (b0'*b1 + b1'*b0)*v7*d7 + (b0'*b2 + b1'*b1 + b2'*b0)*v7*d7.^2 + (b1'*b2 + b2'*b1)*v7*d7.^3 + b2'*b2*v7*d7.^4

