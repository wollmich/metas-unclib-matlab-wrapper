function x = optimizer(func, xStart, p, varargin)
% LinProp Optimizer
%   x = OPTIMIZER(func, xStart, p)
%   x = OPTIMIZER(func, xStart, p, covarianceWeighting)
%   x = OPTIMIZER(func, xStart, p, covarianceWeighting, weights)
%   x = OPTIMIZER(func, xStart, p, covarianceWeighting, weights, bndL, bndU)
%   x = OPTIMIZER(func, xStart, p, covarianceWeighting, weights, bndL, bndU, epsx)
%   x = OPTIMIZER(func, xStart, p, covarianceWeighting, weights, bndL, bndU, epsx, algorithm)

% Michael Wollensack METAS - 29.09.2017

% x Start Values
xStart = double(xStart);
% Parameters
p = LinProp(p);
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


global LinPropOptimizerObjectiveFunction;
LinPropOptimizerObjectiveFunction = func;

xStart2 = Convert2NET(LinProp(xStart));
p2 = Convert2NET(LinProp(p));

d = NET.createGeneric('Metas.UncLib.Optimization.ObjectiveFunction', {'Metas.UncLib.LinProp.UncNumber'}, @ObjectiveFunction2);
x2 = NET.invokeGenericMethod('Metas.UncLib.Optimization.Optimizer', 'Start', {'Metas.UncLib.LinProp.UncNumber'}, d, xStart2, p2, weights, bndL, bndU, epsx, [], [], true, covarianceWeighting, [], 0, algorithm2);

x = Convert2LinProp(x2);
end

function f = ObjectiveFunction2(x, p, debug)

global LinPropOptimizerObjectiveFunction;

x2 = Convert2LinProp(x);
p2 = Convert2LinProp(p);
f2 = LinPropOptimizerObjectiveFunction(x2(:), p2(:));
f = Convert2NET(f2(:));

end

function f = Convert2NET(x)

n = numel(x);
if n == 1
    f = NET.createArray('Metas.UncLib.LinProp.UncNumber', n);
    t = x(1);
    f(1) = t.NetObject;
else
    t = get_net_object(x);
    f = t.data;
end
end

function f = Convert2LinProp(x)

n = x.Length;
if n == 1
    f = LinProp(x(1));
else
    t = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
    t.Init1dData(x, false);
    t.Reshape(int32([n 1]));
    f = LinProp(t);
end
end
