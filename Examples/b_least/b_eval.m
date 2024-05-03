% B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 27.10.2023

function px = b_eval(py, b, b_func)
% b_eval evaluates the fit function b_func with the coefficients b at the
% points py.

px = b_func(py, b);
end
