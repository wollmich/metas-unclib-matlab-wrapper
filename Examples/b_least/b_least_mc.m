% B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 14.11.2023

function [b, fweighted, res] = b_least_mc(py, px, b_func, nsim)
% b_least_mc fits the coefficients b of the fit function b_func using Monte
% Carlo and the points py and px.

% make sure that py and px ae column vectors
py = py(:);
px = px(:);
n = length(py);
py_val = double(py);
px_val = double(px);
bStart = b_start(py_val, px_val, b_func);
% optimization
xOptStart = [py_val; bStart];
xOptScale = xOptStart;
xOptScale(xOptStart == 0) = 1;
xOptStart2 = xOptStart./xOptScale;
pOpt = [py; px];

MCPropGlobalN(nsim);
pOptMC = LinProp2MCProp(pOpt);

% covariance weighting
% CV = R'*R
% W*W' = inv(CV) = inv(R'*R) = inv(R)*inv(R')
% W = inv(R)
W = inv(chol(get_covariance(pOpt)));
objective = @(x2, p2) b_objective(x2, p2, xOptScale, b_func, W);
xOpt2 = optimizer_mc(objective, xOptStart2, pOptMC, false, [], [], [], 1e-9);
xOpt = xOpt2.*xOptScale;
b = xOpt(n+1:end);
% residual
fweighted = objective(xOpt2, pOptMC);
res = fweighted*fweighted';
end

% Start values
function bStart = b_start(y, x, b_func)
if isequal(b_func, @b_linear_func)
    bStart = flip(polyfit(y, x, 1));
    bStart = bStart(:);
elseif isequal(b_func, @b_second_order_ploy)
    bStart = flip(polyfit(y, x, 2));
    bStart = bStart(:);
elseif isequal(b_func, @b_thrid_order_ploy)
    bStart = flip(polyfit(y, x, 3));
    bStart = bStart(:);
elseif isequal(b_func, @b_power_func)
    b1_b2 = flip(polyfit(y, x, 1));
    bStart = [b1_b2(:); 1];
elseif isequal(b_func, @b_exp_func)
    % x = b1 + b2.*exp(b3.*y)
    % x = b1 + b2.*(1 + b3.*y + b3.^2./2.*y.^2 + ...)
    % x = b1 + b2 + b2.*b3.*y + b2.*b3.^2./2.*y.^2 + ...
    % x = c1 + c2.*y + c3.*y.^2 + ...
    c = flip(polyfit(y, x, 2));
    b3 = 2.*c(3)./c(2);
    b2 = c(2)./b3;
    b1 = c(1) - b2;
    bStart = [b1; b2; b3];
else
    error('Unknown fit function')
end
end

% Objective function
% distance from a point to a function
function f = b_objective(x, p, x_scale, b_func, W)
% make sure that x and p are column vectors
x = x(:).*x_scale;
p = p(:);
n = int32(length(p)./2);
py = p(1:n);
px = p(n+1:end);
fy = x(1:n);
b = x(n+1:end);
fx = b_func(fy, b);
dy = py - fy;
dx = px - fx;
f = [dy; dx]'*W;
end
