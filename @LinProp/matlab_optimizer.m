function x = matlab_optimizer(func, xStart, p)
% Matlab Optimizer
% Johannes Hoffmann METAS - 25.02.2013
% Michael Wollensack METAS - 05.03.2013

xStart = double(xStart);
p = LinProp(p);

%fun=@(x) func(x,get_value(p));
fun=@(x) wfunc(func,x,get_value(p),p);

lb=-Inf*ones(1,length(xStart));
ub=Inf*ones(1,length(xStart));
options=optimset('Display','off','MaxFunEvals',300000,'MaxIter',200,'MaxPCGIter',1000,'LargeScale','on','Tolfun',1e-12,'TolX',1e-12,'TolPCG',1e-13,'DiffMaxChange',1e-5,'DiffMinChange',1e-15);
result=lsqnonlin(fun,xStart,lb,ub,options);

x_unc=LinProp(result,diag(ones(length(result),1)));
Jfx=get_jacobi2(fun(x_unc),x_unc);
p_unc=LinProp(get_value(p),diag(ones(length(p),1)));
%fun=@(x) func(x,p_unc);
fun=@(x) wfunc(func,x,p_unc,p);
Jfp=get_jacobi2(fun(result),p_unc);
Jxf=inv(Jfx'*Jfx)*Jfx';
Jxp=-Jxf*Jfp;
nf = length(result);
x = LinProp(zeros(nf,1));
x_net = get_net_object(x);
for i=1:nf
   %x(i) = LinProp(result(i), p, Jxp(i,:), 'system');
   xi = LinProp(result(i), p, Jxp(i,:), 'system');
   xi_net = get_net_object(xi);
   x_net.data(i) = xi_net;
end
end

function f2 = wfunc(func, x2, p2, p)

f = func(x2, p2);
x = LinProp(double(x2));
f_cov = get_covariance(func(x, p));
[V,D] = eig(f_cov);
dd = zeros(size(D));
for i = 1:size(D,1)
    if (D(i,i) > 1e-15)
        dd(i,i) = 1/sqrt(D(i,i));
    end
end
inv_f_jac = V * dd;
f2 = f * inv_f_jac;
end
