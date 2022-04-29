function X = idft(A, varargin)
% LinProp idft
%   X = idft(A) or idft(A, realTimeDomainData) computes the inverse discret fourier transform.

% Michael Wollensack METAS - 29.04.2022

    numlib = NET.createGeneric('Metas.UncLib.LinProp.Ndims.ComplexUncNumLib', {'Metas.UncLib.LinProp.UncNumber'});  
    if (length(varargin) == 0)
        realTimeDomainData = false;
    else
        realTimeDomainData = varargin{1};
    end
    A = complex(A);
    s = size(A);
    am = LinProp.Convert2UncArray(A);
    xm = numlib.Idft(am, realTimeDomainData);
    X = LinProp.Convert2LinProp(xm);
    if (realTimeDomainData)
        X = real(X);
    end
    X = reshape(X, s);
end
