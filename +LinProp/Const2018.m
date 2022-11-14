classdef Const2018
% Physical Constants CODATA 2018
    
    properties (Constant)
        deltavCs =  LinProp.Const2018.data.deltavCs;    % Hyperfine transition frequency of Cs-133 / Hz
        c0 =        LinProp.Const2018.data.c0;          % Speed of light in vacuum / (m/s)
        h =         LinProp.Const2018.data.h;           % Planck constant / Js
        e =         LinProp.Const2018.data.e;           % Elementary charge / C
        k =         LinProp.Const2018.data.k;           % Boltzmann constant / (J/K)
        Na =        LinProp.Const2018.data.Na;          % Avogadro constant / (1/mol)
        Kcd =       LinProp.Const2018.data.Kcd;         % Luminous efficacy / (lm/W)
        Kj =        LinProp.Const2018.data.Kj;          % Josephson constant / (Hz/V)
        Rk =        LinProp.Const2018.data.Rk;          % von Klitzing constant / Ohm
        F =         LinProp.Const2018.data.F;           % Faraday constant / (C/mol)
        R =         LinProp.Const2018.data.R;           % Molar gas constant / (J/(mol*K))
        eV =        LinProp.Const2018.data.eV;          % Electron volt / J
        G =         LinProp.Const2018.data.G;           % Newtonian constant of gravitation / (m^3/(kg*s^2))
        alpha =     LinProp.Const2018.data.alpha;       % Fine-structure constant
        mu0 =       LinProp.Const2018.data.mu0;         % Vacuum magnetic permeability / (Vs/Am)
        ep0 =       LinProp.Const2018.data.ep0;         % Vacuum electric permittivity / (As/Vm)
        Ryd =       LinProp.Const2018.data.Ryd;         % Rydberg constant / (1/m)
        me =        LinProp.Const2018.data.me;          % Electron mass / kg
        are =       LinProp.Const2018.data.are;         % Electron relative atomic mass
        arp =       LinProp.Const2018.data.arp;         % Proton relative atomic mass
        mpsme =     LinProp.Const2018.data.mpsme;       % Proton-electron mass ratio
        mp =        LinProp.Const2018.data.mp;          % Proton mass / kg
        u =         LinProp.Const2018.data.u;           % Atomic mass constant / kg
        Mu =        LinProp.Const2018.data.Mu;          % Molar mass constant / (kg/mol)
    end
    
    properties (Constant, Hidden, Access = private)
        data = loadData(); % Internal field used for loading data - Required to ensure UncPropLoadNETAssemblies('LinProp'); is called
    end
    
end

function data = loadData()
    UncPropLoadNETAssemblies('LinProp');
    
    data = struct();
    
    data = addFields(data, Metas.UncLib.Core.Const2018);
    data = addFields(data, NET.createGeneric('Metas.UncLib.Core.Const2018', {'Metas.UncLib.LinProp.UncNumber'}), @(x) LinProp(x));

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
