% B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 27.10.2023

function x = b_thrid_order_ploy(y, b)
x = b(1) + b(2).*y + b(3).*y.^2 + b(4).*y.^3;
end
