function c = get_expanded_covariance(x, p)
% GET_EXPANDED_FACTOR returns the expanded covariance matrix

% Michael Wollensack METAS - 15.03.2023

x = MCProp(x);
cv = get_covariance(x);
p = get_coverage_factor(x, cv, p);
c = p.*p.*cv;
end
