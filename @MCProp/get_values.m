function values = get_values(x)
% GET_VALUES returns the Monte Carlo random numbers 

% Michael Wollensack METAS - 15.03.2023

x = MCProp(x);
l = ToUncList(x);
n1 = MCPropGlobalN;
n2 = l.data.Length;
values = zeros(n1, n2);
for i2 = 1:n2
    xi2n = l.data(i2);
    values(:,i2) = double(xi2n.values);
end

end
