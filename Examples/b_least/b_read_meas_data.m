% B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 07.12.2023

function pym = b_read_meas_data(filepath)
% b_read_meas_data reads measurement data from tabular separated text file
% where the first column are the y values and the second column are the
% standard uncertainties of y.

d = dlmread(filepath);
n = size(d, 1);
pym = LinProp(zeros(n, 1));
for i = 1:n
    pym(i) = LinProp(d(i, 1), d(i, 2), ['Measurement No. ' num2str(i) ' y']);
end
end
