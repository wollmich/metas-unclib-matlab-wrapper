function dY = Example_ODE_Function(Y, x, P, dY)
% Example ODE Function
% Michael Wollensack METAS - 14.10.2025

% ODE Function f(Y, x, P, dY)

dY(1) = P(1).*Y(1);

end
