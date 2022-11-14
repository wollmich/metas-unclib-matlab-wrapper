classdef Const2014
% Physical Constants CODATA 2014
    
    properties (Constant)
        deltavCs =  MCProp.Const2014.data.deltavCs;    % Hyperfine transition frequency of Cs-133 / Hz
        c0 =        MCProp.Const2014.data.c0;          % Speed of light in vacuum / (m/s)
        mu0 =       MCProp.Const2014.data.mu0;         % Vacuum magnetic permeability / (Vs/Am)
        ep0 =       MCProp.Const2014.data.ep0;         % Vacuum electric permittivity / (As/Vm)
        Kcd =       MCProp.Const2014.data.Kcd;         % Luminous efficacy / (lm/W)
        Mu =        MCProp.Const2014.data.Mu;          % Molar mass constant / (kg/mol)
        G =         MCProp.Const2014.data.G;           % Newtonian constant of gravitation / (m^3/(kg*s^2))
        alpha =     MCProp.Const2014.data.alpha;       % Fine-structure constant
        Ryd =       MCProp.Const2014.data.Ryd;         % Rydberg constant / (1/m)
        mpsme =     MCProp.Const2014.data.mpsme;       % Proton-electron mass ratio
        Na =        MCProp.Const2014.data.Na;          % Avogadro constant / (1/mol)
        Kj =        MCProp.Const2014.data.Kj;          % Josephson constant / (Hz/V)
        k =         MCProp.Const2014.data.k;           % Boltzmann constant / (J/K)
        Rk =        MCProp.Const2014.data.Rk;          % von Klitzing constant / Ohm
        e =         MCProp.Const2014.data.e;           % Elementary charge / C
        h =         MCProp.Const2014.data.h;           % Planck constant / Js
        me =        MCProp.Const2014.data.me;          % Electron mass / kg
        mp =        MCProp.Const2014.data.mp;          % Proton mass / kg
        u =         MCProp.Const2014.data.u;           % Atomic mass constant / kg
        F =         MCProp.Const2014.data.F;           % Faraday constant / (C/mol)
        R =         MCProp.Const2014.data.R;           % Molar gas constant / (J/(mol*K))
        eV =        MCProp.Const2014.data.eV;          % Electron volt / J
    end
    
    properties (Constant, Access = private)
        data = loadData(); % Internal field used for loading data - Required to ensure UncPropLoadNETAssemblies('MCProp'); is called
    end
    
end

function data = loadData()
    UncPropLoadNETAssemblies('MCProp');
    
    data = struct();
    
    data = addFields(data, Metas.UncLib.Core.Const2014);
    data = addFields(data, NET.createGeneric('Metas.UncLib.Core.Const2014', {'Metas.UncLib.MCProp.UncNumber'}), @(x) MCProp(x));

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


