function Y = Example_Numerical_Function(X)
% Example Numerical Function
% Michael Wollensack METAS - 28.02.2013

% Numerical Function Y = f(X)

disp(X(:).');
 
a = X(1);
b = X(2);
c = sqrt(a.*a + b.*b);
d = a + b + c;
Y = [c d];

end
