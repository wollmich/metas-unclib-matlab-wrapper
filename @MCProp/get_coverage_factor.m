function k = get_coverage_factor(x, cv, p)
% GET_COVERAGE_FACTOR returns the coverage factor

% Michael Wollensack METAS - 15.03.2023

v = get_values(x);
n1 = MCPropGlobalN;
e = zeros(n1, 1);
icv = inv(cv);
for i1 = 1:n1
    e(i1) = sqrt(v(i1,:)*icv*v(i1,:)');
end
f = sort(e);
k = f(round(p.*n1));
end
