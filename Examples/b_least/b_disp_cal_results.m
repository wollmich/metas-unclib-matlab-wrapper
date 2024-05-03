% B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 07.12.2023

function b_disp_cal_results(b, fweighted, res)

cov = get_covariance(b);
fprintf('\ncoefficients b ± u(b)\n')
for i = 1:length(b)
    fprintf('b(%d) = %s\n', i, b(i))
end
for i = 1:length(b)
    for j = i+1:length(b)
        fprintf('covariance between b(%d) and b(%d) = %0.4e\n', j, i, cov(j, i))
    end
end
fprintf('residual  = %0.4f\n', res)
fprintf('maximum absolute value of weighted deviations = %0.4f\n', max(abs(fweighted)))
end
