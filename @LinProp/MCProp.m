function xMC = MCProp(x)
% Converts objects to MCProp
%
%   xMC = MCProp(x) returns MCProp objects where
%     x are the input LinProp objects.
%
%   Example of usage:
%     xMC = MCProp(x);
%     yMC = f(xMC);
%     y = MCProp2LinProp(yMC, xMC, x);
%
%   The expected values of y are the same as the expected values of yMC.
%   The covariance of y is the same as the covariance of yMC.

% Michael Wollensack METAS - 22.02.2023

x = LinProp(x);
xValues = get_value(x);
xCovariance = get_covariance(x);

if numel(xCovariance) == 1
    xMC = MCProp(xValues, sqrt(xCovariance));
else
    xMC = MCProp(xValues, xCovariance);
end

end
