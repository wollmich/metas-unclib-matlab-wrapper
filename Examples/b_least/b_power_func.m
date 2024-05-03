% B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 07.12.2023

function x = b_power_func(y, b)
x = b(1) + b(2).*y.^(1 + b(3));
end
