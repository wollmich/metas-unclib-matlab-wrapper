function x = optimizer_mc(func, xStart, p, varargin)
% MCProp Optimizer
%   x = OPTIMIZER_MC(func, xStart, p)
%   x = OPTIMIZER_MC(func, xStart, p, covarianceWeighting)
%   x = OPTIMIZER_MC(func, xStart, p, covarianceWeighting, weights)
%   x = OPTIMIZER_MC(func, xStart, p, covarianceWeighting, weights, bndL, bndU)
%   x = OPTIMIZER_MC(func, xStart, p, covarianceWeighting, weights, bndL, bndU, epsx)
%   x = OPTIMIZER_MC(func, xStart, p, covarianceWeighting, weights, bndL, bndU, epsx, algorithm)

% Michael Wollensack METAS - 14.11.2023

% x Start Values
xStart = double(xStart);
% Parameters
p = MCProp(p);
% Covariance Weighting
covarianceWeighting = true;
if (length(varargin) > 0) covarianceWeighting = logical(varargin{1}); end
% Weights
weights = [];
if (length(varargin) > 1) weights = double(varargin{2}); end
% Lower Bounds
bndL = [];
if (length(varargin) > 2) bndL = double(varargin{3}); end
% Upper Bounds
bndU = [];
if (length(varargin) > 3) bndU = double(varargin{4}); end
% Eps x
epsx = 0;
if (length(varargin) > 4) epsx = double(varargin{5}); end
% Algorithm
algorithm = 'LevenbergMarquardt'; % or 'TrustRegion'
if (length(varargin) > 5) algorithm = char(varargin{6}); end
algorithm2 = eval(['Metas.UncLib.Optimization.Algorithm.' algorithm]);

if (covarianceWeighting)
    error('Covariance weighting not supported');
end

pval = get_values(p(:)');
n1 = size(pval, 1);
n2 = length(xStart);
xval = zeros(n1, n2);
for i1 = 1:n1
    pval_i1 = LinProp(pval(i1,:));
    x_i1 = optimizer(func, xStart, pval_i1, false, weights, bndL, bndU, epsx, algorithm);
    xval(i1,:) = double(x_i1(:))';
end
x = MCProp(zeros(n2, 1));
for i2 = 1:n2
    x_i2_n = get_net_object(x(i2));
    x_i2_n.fcn_value = mean(xval(:,i2));
    x_i2_n.values = xval(:,i2);
    x(i2) = MCProp(x_i2_n);
end
