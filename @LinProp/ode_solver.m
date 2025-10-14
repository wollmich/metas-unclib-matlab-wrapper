function ytbl = ode_solver(y, x, p, eps, func)
% LinProp ODE Solver
% Michael Wollensack METAS - 14.10.2025

y = LinProp(y);
x = LinProp(x);
p = LinProp(p);

global LinPropOdeFunction;
LinPropOdeFunction = func;

y2 = Convert2NET(y);
x2 = Convert2NET(x);
p2 = Convert2NET(p);

d = Metas.UncLib.LinProp.OdeFunction2(@Temp);
ytbl2 = NET.invokeGenericMethod('Metas.UncLib.LinProp.UncOdeSolver', 'Solve', {'Metas.UncLib.LinProp.UncNumber'}, y2, x2, p2, eps, d);

ytbl = Convert2LinProp(ytbl2);
end

function y = Temp(y, x, p, dy)

global LinPropOdeFunction;

y2 = double(y);
x2 = double(x);
p2 = double(p);
dy2 = double(dy);
dy2 = LinPropOdeFunction(y2, x2, p2, dy2);
for i = 1:length(dy2)
    dy(i) = dy2(i);
end
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

t = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
t.Init2dData(x);
f = LinProp(t);
end
