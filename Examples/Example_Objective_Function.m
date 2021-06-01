function F = Example_Objective_Function(X, P)
% Example Objective Function
% Michael Wollensack METAS - 28.02.2013

% Objective function F = f(X, P)
% X : variable optimization parameters
% P : constant optimization parameters

disp(double(X(:)).');

f1 = P(1)*X(1)*X(1) + P(2)*X(1) + P(3);
f2 = exp(X(2)) - P(4);
F = [f1 f2];

end
