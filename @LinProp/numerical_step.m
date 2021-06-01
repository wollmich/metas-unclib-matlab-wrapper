function y = numerical_step(func, x, dx)
% LinProp Monte Carlo Step
% Michael Wollensack METAS - 22.02.2017

x = LinProp(x);
dx = double(dx);

global LinPropNumericalStepFunction;
LinPropNumericalStepFunction = func;

x2 = Convert2NET(x);

%n = Metas.UncLib.Core.Parallel.ThreadCount;
%NET.setStaticProperty('Metas.UncLib.Core.Parallel.ThreadCount', 1)

d = Metas.UncLib.LinProp.NumericalFunctionDelegate(@Temp);
y2 = Metas.UncLib.LinProp.UncNumerical.NumericalStep(d, x2, dx, true);

%NET.setStaticProperty('Metas.UncLib.Core.Parallel.ThreadCount', n)

y = Convert2LinProp(y2);
end

function y = Temp(x)

global LinPropNumericalStepFunction;

x = double(x);
y = LinPropNumericalStepFunction(x);

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
t = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
t.Init1dData(x, false);
t.Reshape(int32([n 1]));
f = LinProp(t);
end
