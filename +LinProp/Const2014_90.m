classdef Const2014_90
% Physical Constants CODATA 2014 for Conventional Electrical Units 90
    
    properties (Constant)
        deltavCs =  LinProp.Const2014.data.deltavCs;    % Hyperfine transition frequency of Cs-133 / Hz
        c0 =        LinProp.Const2014.data.c0;          % Speed of light in vacuum / (m/s)
        mu0 =       LinProp.Const2014.data.mu0;         % Vacuum magnetic permeability / (Vs/Am)
        ep0 =       LinProp.Const2014.data.ep0;         % Vacuum electric permittivity / (As/Vm)
        Kcd =       LinProp.Const2014.data.Kcd;         % Luminous efficacy / (lm/W)
        Mu =        LinProp.Const2014.data.Mu;          % Molar mass constant / (kg/mol)
        Kj =        LinProp.Const2014.data.Kj;          % Conventional value of Josephson constant / (Hz/V)
        Rk =        LinProp.Const2014.data.Rk;          % Conventional value of von Klitzing constant / Ohm
        e =         LinProp.Const2014.data.e;           % Elementary charge / C
        h =         LinProp.Const2014.data.h;           % Planck constant / Js
        Na =        LinProp.Const2014.data.Na;          % Avogadro constant / (1/mol)
        F =         LinProp.Const2014.data.F;           % Faraday constant / (C/mol)
        k =         LinProp.Const2014.data.k;           % Boltzmann constant / (J/K)
    end
    
    properties (Constant, Access = private)
        data = loadData(); % Internal field used for loading data - Required to ensure UncPropLoadNETAssemblies('LinProp'); is called
    end
    
end

function data = loadData()
    UncPropLoadNETAssemblies('LinProp');
    
    data = struct();
    
    data = addFields(data, Metas.UncLib.Core.Const2014);
    data = addFields(data, Metas.UncLib.Core.Const2014_90);
    data = addFields(data, NET.createGeneric('Metas.UncLib.Core.Const2014_90', {'Metas.UncLib.LinProp.UncNumber'}), @(x) LinProp(x));

end

function s = addFields(s, s2, fH)
% Adds all fields of s2 to s. If defined, the function handle fH is
% applied to each field of s2.
    if nargin < 3
        fH = @(x) x;
    end

    fields = fieldnames(s2);
    for ii = 1:numel(fields)
        s.(fields{ii}) = fH(s2.(fields{ii}));
    end
end
