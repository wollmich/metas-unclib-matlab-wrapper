function [V,D] = eig(A0, varargin)
% LinProp Eigenvalues and eigenvectors
%   [V,D] = EIG(A0)
%
%   [V,D] = EIG(A0, A1) linear Eigenvalue problem A0*V + A1*V*D = 0
%
%   [V,D] = EIG(A0, A1, A2)quadratic Eigenvalue problem A0*V + A1*V*D + A2*V*D.^2 = 0
%
%   [V,D] = EIG(A0, A1, A2, A3) cubic Eigenvalue problem A0*V + A1*V*D + A2*V*D.^2 + A3*V*D.^3 = 0
%
%   [V,D] = EIG(A0, A1, A2, A3, A4) quartic Eigenvalue problem A0*V + A1*V*D + A2*V*D.^2 + A3*V*D.^3 + A4*V*D.^4 = 0
%
%   [V,D] = EIG(A0, A1, A2, ..., An) non-linear Eigenvalue problem A0*V + A1*V*D + A2*V*D.^2 + ... + An*V*D.^n = 0

% Michael Wollensack METAS - 20.08.2018

A0 = LinProp(A0);
isComplex = A0.IsComplex;
for i = 1:length(varargin)
    varargin{i} = LinProp(varargin{i});
    isComplex = isComplex | varargin{i}.IsComplex;
end

n = size(A0);
n1 = n(1);
n2 = n(2);

if (length(varargin) == 0)
    if (n1 ~= n2)
        throw('Matrix must be square');
    end
    if (~isComplex)
        symmetric = true;
        A0n = Convert2UncArray(A0);
        for i1 = 1:n(1)
            for i2 = i1 + 1:n(2)
                if (symmetric)
                    diff = LinProp(A0n.GetItem2d(i1-1, i2-1).Subtract(A0n.GetItem2d(i2-1, i1-1)));
                    symmetric = symmetric & get_value(diff) == 0 & get_stdunc(diff) == 0;
                end
            end
        end
        isComplex = ~symmetric;
    end
else
    isComplex = true;
end

if (isComplex)
    A0 = complex(A0);
    for i = 1:length(varargin)
        varargin{i} = complex(varargin{i});
    end
end

linalg = UncLinAlg(isComplex);
if (length(varargin) == 0)
    A0n = Convert2UncArray(A0);
    if (isComplex)
        [Vn,Dn] = linalg.NonsymmetricEig(A0n);
    else
        [Vn,Dn] = linalg.SymmetricEig(A0n);
    end
else
    An = InitArray(isComplex, length(varargin) + 1);
    An(1) = Convert2UncArray(A0);
    for i = 1:length(varargin)
        An(i + 1) = Convert2UncArray(varargin{i});
    end
    [Vn,Dn] = linalg.NonLinearEig(An);
end

Vn2 = SubMatrix(Vn, isComplex, n2);
Vn3 = RemoveZeroImagPart(Vn2, isComplex);
Dn2 = DiagMatrix(Dn, isComplex);
Dn3 = RemoveZeroImagPart(Dn2, isComplex);
V = Convert2LinProp(Vn3);
D = Convert2LinProp(Dn3);
end

function m = InitArray(isComplex, n)
    if isComplex
        c = NET.GenericClass('Metas.UncLib.Core.Ndims.ComplexNArray', 'Metas.UncLib.LinProp.UncNumber');
    else
        c = NET.GenericClass('Metas.UncLib.Core.Ndims.RealNArray', 'Metas.UncLib.LinProp.UncNumber');
    end
    m = NET.createArray(c, n);
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

function m = DiagMatrix(d, isComplex)
    if isComplex
        m = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
    else
        m = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
    end
    n = d.numel;
    m.Zeros2d(n, n);
    for i=1:n
        m.SetItem2d(i-1,i-1,d.GetItem1d(i-1));
    end
end

function m = SubMatrix(x, isComplex, nlines)
    n = x.size;
    if n(1) > nlines
        if isComplex
            m = NET.createGeneric('Metas.UncLib.Core.Ndims.ComplexNArray', {'Metas.UncLib.LinProp.UncNumber'});
        else
            m = NET.createGeneric('Metas.UncLib.Core.Ndims.RealNArray', {'Metas.UncLib.LinProp.UncNumber'});
        end
        m.Zeros2d(nlines, n(2));
        for i1=1:nlines
            for i2=1:n(2)
                m.SetItem2d(i1-1,i2-1,x.GetItem2d(i1-1,i2-1));
            end
        end
    else
        m = x;
    end
end

function m = RemoveZeroImagPart(x, isComplex)
    if isComplex
        n = x.numel;
        zeroImag = true;
        for i = 1:n
            if zeroImag
                xi = x.GetItem1d(i-1);
                zeroImag = zeroImag & xi.imag.IsZero;
            end
        end
        if zeroImag
            m = x.Real();
        else
        m = x;
        end
    else
    m = x;
    end
end

function l = UncLinAlg(complex)
    if complex
        l = NET.createGeneric('Metas.UncLib.LinProp.Ndims.ComplexUncLinAlg', {'Metas.UncLib.LinProp.UncNumber'});
    else
        l = NET.createGeneric('Metas.UncLib.LinProp.Ndims.RealUncLinAlg', {'Metas.UncLib.LinProp.UncNumber'});
    end
end
