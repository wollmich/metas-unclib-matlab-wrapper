classdef Const2014_90
% Physical Constants CODATA 2014 for Conventional Electrical Units 90
    
    properties (Constant)
        deltavCs =  DistProp.Const2014_90.data.deltavCs;    % Hyperfine transition frequency of Cs-133 / Hz
        c0 =        DistProp.Const2014_90.data.c0;          % Speed of light in vacuum / (m/s)
        mu0 =       DistProp.Const2014_90.data.mu0;         % Vacuum magnetic permeability / (Vs/Am)
        ep0 =       DistProp.Const2014_90.data.ep0;         % Vacuum electric permittivity / (As/Vm)
        Kcd =       DistProp.Const2014_90.data.Kcd;         % Luminous efficacy / (lm/W)
        Mu =        DistProp.Const2014_90.data.Mu;          % Molar mass constant / (kg/mol)
        Kj =        DistProp.Const2014_90.data.Kj;          % Conventional value of Josephson constant / (Hz/V)
        Rk =        DistProp.Const2014_90.data.Rk;          % Conventional value of von Klitzing constant / Ohm
        e =         DistProp.Const2014_90.data.e;           % Elementary charge / C
        h =         DistProp.Const2014_90.data.h;           % Planck constant / Js
        Na =        DistProp.Const2014_90.data.Na;          % Avogadro constant / (1/mol)
        F =         DistProp.Const2014_90.data.F;           % Faraday constant / (C/mol)
        k =         DistProp.Const2014_90.data.k;           % Boltzmann constant / (J/K)
    end
    
    properties (Constant, Access = private)
        data = loadData(); % Internal field used for loading data - Required to ensure UncPropLoadNETAssemblies('DistProp'); is called
    end
    
end

function data = loadData()
    UncPropLoadNETAssemblies('DistProp');
    
    data = struct();
    
    data = addFields(data, Metas.UncLib.Core.Const2014);
    data = addFields(data, Metas.UncLib.Core.Const2014_90);
    data = addFields(data, NET.createGeneric('Metas.UncLib.Core.Const2014_90', {'Metas.UncLib.DistProp.UncNumber'}), @(x) DistProp(x));

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
