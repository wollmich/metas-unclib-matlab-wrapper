% B_LEAST ISO 6143:2001 
% Michael Wollensack METAS - 25.10.2023 - 01.11.2023

function h = b_plot(py, px, pym, pxm, b, b_func)
% b_plot plots the reference points py and pxm, the fit function b_func
% uisng the coefficients b and the measurement points pym and pxm.

k = 2;
h = figure();
% fit
ymin = min([double(py(:)); double(pym(:))]);
ymax = max([double(py(:)); double(pym(:))]);
fy = linspace(ymin, ymax, 100);
fx = b_eval(fy, b, b_func);
% fit uncertainty region
Ufx = k.*get_stdunc(fx);
fx_p = double(fx) + Ufx;
fx_n = double(fx) - Ufx;
fill([fx_p(1:end) fx_n(end:-1:1)], [fy(1:end) fy(end:-1:1)], 'b', 'FaceAlpha', 0.5, 'LineStyle','none');
hold on;
% uncertainty region of reference points
for i=1:length(py)
    plot_ellipse(px(i), py(i), 'r')
end
% uncertainty region of measurement points
for i=1:length(pym)
    plot_ellipse(pxm(i), pym(i), 'k')
end
% reference points
Upy = k.*get_stdunc(py);
Upx = k.*get_stdunc(px);
hp1 = errorbar(double(px), double(py), Upy, Upy, Upx, Upx, 'r+', 'LineWidth', 1, 'DisplayName', 'Reference points');
% fit function
hp2 = plot(double(fx), fy, 'b-', 'LineWidth', 1, 'DisplayName', 'Fit x=f(y)');
% measurement points
Upym = k.*get_stdunc(pym);
Upxm = k.*get_stdunc(pxm);
hp3 = errorbar(double(pxm), double(pym), Upym, Upym, Upxm, Upxm, 'k+', 'LineWidth', 1, 'DisplayName', 'Measurement points');
% title, labels and legend
hold off;
xlabel('Assigned value x');
ylabel('Instrument response y');legend([hp1 hp2 hp3], 'Location', 'best');
grid on;
width = 800; 
height = 600;
Pix_SS = get(0,'screensize');
fig=gcf;
fig.Position = [(Pix_SS(3)-width)/2 (Pix_SS(4)-height)/2 width height];
end

function plot_ellipse(pxi, pyi, c)
k = 2.45;
cv = get_covariance([pxi pyi]);
[V, D] = eig(k.*cv);
t = linspace(0, 2*pi, 100);
e = V*sqrt(D)*[cos(t); sin(t)];
fill(e(1,:) + double(pxi), e(2,:) + double(pyi), c, 'FaceAlpha', 0.5, 'LineStyle','none');
end
