% B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 27.10.2023

function b_disp_meas_results(pxm, pym)

fprintf('\nprediction x ± u(x) based on y ± u(y)\n')
for i = 1:length(pxm)
    fprintf('Measurement No. %d x = %s\t y = %s\n', i, pxm(i), pym(i))
end
