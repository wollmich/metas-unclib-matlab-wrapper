function X = dft(A)
% LinProp dft
%   X = dft(A) computes the discret fourier transform.

% Michael Wollensack METAS - 29.04.2022

    numlib = NET.createGeneric('Metas.UncLib.LinProp.Ndims.ComplexUncNumLib', {'Metas.UncLib.LinProp.UncNumber'});
    realTimeDomainData = ~A.IsComplex;
    A = complex(A);
    s = size(A);
    am = LinProp.Convert2UncArray(A);
    xm = numlib.Dft(am, realTimeDomainData);
    X = LinProp.Convert2LinProp(xm);
    X = reshape(X, s);
end
