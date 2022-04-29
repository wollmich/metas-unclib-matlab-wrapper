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
    am = Convert2UncArray(A);
    xm = numlib.Idft(am, realTimeDomainData);
    X = Convert2LinProp(xm);
    if (realTimeDomainData)
        X = real(X);
    end
    X = reshape(X, s);
end

function m = Convert2UncArray(x)
    if x.IsArray
        m = x.NetObject;
    else
        if x.IsComplex
            m = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
        else
            m = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
        end
        m.Init2d(1, 1);
        m.SetItem2d(0, 0, x.NetObject);
    end 
end

function u = Convert2LinProp(x)
    if LinProp.IsArrayNet(x)
        if x.numel == 1
            u = LinProp(x.GetItem2d(0, 0));
        else
            u = LinProp(x);
            if ndims(u) == 1
                u = reshape(u, size(u));
            end
        end
    else
        u = LinProp(x);
    end
end
