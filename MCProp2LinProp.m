function y = MCProp2LinProp(yMC, xMC, x)
% MCProp2LinProp Converts MCProp objects back to LinProp objects
%
%   y = MCProp2LinProp(yMC, xMC, x) returns LinProp objects where
%     yMC are the output MCProp objects,
%     xMC are the input MCProp objects and
%     x are the input LinProp objects.
%
%   Example of usage:
%     xMC = LinProp2MCProp(x);
%     yMC = f(xMC);
%     y = MCProp2LinProp(yMC, xMC, x);
%
%   The expected values of y are the same as the expected values of yMC.
%   The covariance of y is the same as the covariance of yMC.

% Michael Wollensack METAS - 02.03.2023

yMC = MCProp(yMC);
xMC = MCProp(xMC);
x = LinProp(x);

if isreal(yMC)
    yMC2 = yMC;
else
    yMC2 = [real(yMC(:)) imag(yMC(:))]';
end
if isreal(xMC) && isreal(x)
    xMC2 = xMC;
    x2 = x;
else
    xMC2 = [real(xMC(:)) imag(xMC(:))]';
    x2 = [real(x(:)) imag(x(:))]';
end
y2 = MCProp2LinPropSub(yMC2, xMC2, x2);
if isreal(yMC)
    y = y2;
else
    y2 = reshape(y2, size(yMC2));
    y = y2(1,:) + 1i.*y2(2,:);
end
y = reshape(y, size(yMC));
end

function y = MCProp2LinPropSub(yMC, xMC, x)

xMCsub = MCProp(zeros(1,0));
xsub = LinProp(zeros(1,0));
n = numel(x);
for i = 1:n
    if get_stdunc(x(i)) > 0
        xMCsub = [xMCsub xMC(i)];
        xsub = [xsub x(i)];
    end
end

xMC_Covar = get_covariance(xMCsub);
x_Covar = get_covariance(xsub);
xMC_x_Jacobi = chol(xMC_Covar, 'lower') * inv(chol(x_Covar,'lower'));

n = numel(xMCsub);
x2 = LinProp(zeros(n, 1));
for i = 1:n
    x2(i) = LinProp(get_value(xMCsub(i)), xsub, xMC_x_Jacobi(i,:), 'system');
end

yMC_xMC_Jacobi = get_jacobi2(yMC, xMCsub);

n = numel(yMC);
y = LinProp(zeros(n, 1));
for i=1:n
    y(i) = LinProp(get_value(yMC(i)), x2, yMC_xMC_Jacobi(i,:), 'system');
end

% Add non-linear contributions
yMC_Covar = get_covariance(yMC);
yMC_Covar_Lin = yMC_xMC_Jacobi*xMC_Covar*yMC_xMC_Jacobi';
yMC_Covar_NonLin = yMC_Covar - yMC_Covar_Lin;

if numel(yMC_Covar_NonLin) == 1
    l = LinProp(0, sqrt(yMC_Covar_NonLin), 'Non-linear contribution');
else
    l = LinProp(zeros(length(yMC_Covar_NonLin), 1), yMC_Covar_NonLin, 'Non-linear contributions');
end
y = y + l;
end
