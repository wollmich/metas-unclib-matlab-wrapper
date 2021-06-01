function r = roots(C)
% LinProp Roots
%   r = roots(C) computes the roots of the polynomial whose coefficients
%    are the elements of the vector C. If C has N+1 components,
%    the polynomial is C(1)*X^N + ... + C(N)*X + C(N+1).

% Michael Wollensack METAS - 20.08.2018

n = numel(C);
if (n < 2)
    r = LinProp(zeros(0,1));
else
    P = num2cell(C(:));
    [~,D] = eig(P{end:-1:1});
    r = diag(D);
end
end
