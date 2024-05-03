% B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 07.12.2023

function [py, px] = b_read_cal_data(filepath)
% b_read_cal_data reads calibration data from tabular separated text file
% where the first column are the x values, the second column are the
% standard uncertainties of x, the third column are the y values and the
% forth column are the standard uncertainties of y.

d = dlmread(filepath);
n = size(d, 1);
py = LinProp(zeros(n, 1));
px = LinProp(zeros(n, 1));
for i = 1:n
    px(i) = LinProp(d(i, 1), d(i, 2), ['Reference No. ' num2str(i) ' x']);
    py(i) = LinProp(d(i, 3), d(i, 4), ['Reference No. ' num2str(i) ' y']);
end
end
